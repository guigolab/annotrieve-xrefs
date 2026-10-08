"""Read-only SQLite session for gene_corpus.sqlite (thread-local)."""
from __future__ import annotations

import sqlite3
import threading
from dataclasses import dataclass
from pathlib import Path

from helpers.paths import gene_corpus_db_path
from settings import settings

_LOCAL = threading.local()
_REGISTRY_LOCK = threading.Lock()
# thread ident -> connection (for process-wide close_reader)
_REGISTRY: dict[int, sqlite3.Connection] = {}

# Soft heap so sorts/hashes spill to temp files instead of growing unbounded.
_SOFT_HEAP_LIMIT = 512 * 1024 * 1024  # 512 MiB


class GeneCorpusUnavailable(Exception):
    """Raised when the published gene corpus DB is missing or unreadable."""


@dataclass
class _ThreadReader:
    path: Path
    mtime: float
    conn: sqlite3.Connection


def _open_connection(resolved: Path) -> sqlite3.Connection:
    uri = f"file:{resolved}?mode=ro&immutable=1"
    # check_same_thread=True: this connection is only used by its owning thread.
    conn = sqlite3.connect(uri, uri=True, check_same_thread=True)
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA cache_size={settings.READ_CACHE_SIZE}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute(f"PRAGMA soft_heap_limit={_SOFT_HEAP_LIMIT}")
    conn.execute("PRAGMA query_only=ON")
    return conn


def close_reader() -> None:
    """Close all thread-local corpus readers (tests / shutdown)."""
    with _REGISTRY_LOCK:
        for conn in _REGISTRY.values():
            try:
                conn.close()
            except sqlite3.Error:
                pass
        _REGISTRY.clear()
    # Drop this thread's cached handle so the next open_reader rebuilds.
    if getattr(_LOCAL, "reader", None) is not None:
        _LOCAL.reader = None


def open_reader(path: Path | None = None) -> sqlite3.Connection:
    """
    Return a thread-local read-only connection to gene_corpus.sqlite.

    Each Uvicorn/gunicorn worker thread keeps its own connection. When the
    published DB path or mtime changes, only **this** thread's connection is
    closed and reopened — other threads refresh on their next call.

    Callers must not close the returned connection; use ``close_reader``.
    """
    db = path if path is not None else gene_corpus_db_path()
    if db is None or not db.is_file():
        raise GeneCorpusUnavailable("Gene corpus is not published")
    resolved = db.resolve()
    mtime = resolved.stat().st_mtime

    state: _ThreadReader | None = getattr(_LOCAL, "reader", None)
    if (
        state is not None
        and state.path == resolved
        and state.mtime == mtime
    ):
        return state.conn

    tid = threading.get_ident()
    if state is not None:
        try:
            state.conn.close()
        except sqlite3.Error:
            pass
        with _REGISTRY_LOCK:
            _REGISTRY.pop(tid, None)

    conn = _open_connection(resolved)
    _LOCAL.reader = _ThreadReader(path=resolved, mtime=mtime, conn=conn)
    with _REGISTRY_LOCK:
        _REGISTRY[tid] = conn
    return conn
