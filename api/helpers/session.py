"""Read-only SQLite session for gene_corpus.sqlite."""
from __future__ import annotations

import sqlite3
import threading
from pathlib import Path

from helpers.paths import gene_corpus_db_path
from settings import settings

_CONN: sqlite3.Connection | None = None
_DB_PATH: Path | None = None
_MTIME: float | None = None
_LOCK = threading.Lock()

# Soft heap so sorts/hashes spill to temp files instead of growing unbounded.
_SOFT_HEAP_LIMIT = 512 * 1024 * 1024  # 512 MiB


class GeneCorpusUnavailable(Exception):
    """Raised when the published gene corpus DB is missing or unreadable."""


def close_reader() -> None:
    """Close the process-wide reader (tests / shutdown / reopen)."""
    global _CONN, _DB_PATH, _MTIME
    with _LOCK:
        if _CONN is not None:
            _CONN.close()
        _CONN = None
        _DB_PATH = None
        _MTIME = None


def open_reader(path: Path | None = None) -> sqlite3.Connection:
    """
    Return a process-wide read-only connection to gene_corpus.sqlite.

    Reopens when the DB path or mtime changes (atomic publish replace).
    Callers must not close the returned connection; use ``close_reader``.
    """
    global _CONN, _DB_PATH, _MTIME
    db = path if path is not None else gene_corpus_db_path()
    if db is None or not db.is_file():
        raise GeneCorpusUnavailable("Gene corpus is not published")
    resolved = db.resolve()
    mtime = resolved.stat().st_mtime
    with _LOCK:
        if _CONN is not None and _DB_PATH == resolved and _MTIME == mtime:
            return _CONN

        if _CONN is not None:
            _CONN.close()
            _CONN = None

        uri = f"file:{resolved}?mode=ro&immutable=1"
        conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        conn.execute(f"PRAGMA cache_size={settings.READ_CACHE_SIZE}")
        conn.execute("PRAGMA temp_store=FILE")
        conn.execute("PRAGMA mmap_size=0")
        conn.execute(f"PRAGMA soft_heap_limit={_SOFT_HEAP_LIMIT}")
        conn.execute("PRAGMA query_only=ON")
        _CONN = conn
        _DB_PATH = resolved
        _MTIME = mtime
        return conn
