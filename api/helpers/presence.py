"""Read-only thread-local connection to ``presence.sqlite`` and taxon-scope bitsets."""
from __future__ import annotations

import sqlite3
import threading
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from helpers.errors import http_error
from helpers.postings import iter_bits  # noqa: F401  (re-export for services)
from settings import presence_dir, settings

PRESENCE_FILE = "presence.sqlite"
_LOCAL = threading.local()
_SOFT_HEAP_LIMIT = 256 * 1024 * 1024


@dataclass
class _Reader:
    path: Path
    mtime: float
    conn: sqlite3.Connection
    all_mask: int
    data_version: str


def presence_enabled() -> bool:
    return presence_dir() is not None


def presence_db_path() -> Path | None:
    root = presence_dir()
    return None if root is None else root / "current" / PRESENCE_FILE


def close_presence() -> None:
    state = getattr(_LOCAL, "reader", None)
    if state is not None:
        try:
            state.conn.close()
        except sqlite3.Error:
            pass
    _LOCAL.reader = None


def _open(resolved: Path) -> _Reader:
    conn = sqlite3.connect(f"file:{resolved}?mode=ro&immutable=1", uri=True,
                           check_same_thread=True)
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA cache_size={settings.READ_CACHE_SIZE}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute(f"PRAGMA soft_heap_limit={_SOFT_HEAP_LIMIT}")
    conn.execute("PRAGMA query_only=ON")
    mask = 0
    for (o,) in conn.execute("SELECT ord FROM annotation"):
        mask |= 1 << int(o)
    meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
    return _Reader(resolved, resolved.stat().st_mtime, conn, mask, str(meta.get("data_version", "")))


def presence_reader() -> _Reader:
    """This thread's reader; reopens when ``current`` points at a new version."""
    db = presence_db_path()
    if db is None or not db.is_file():
        raise http_error(503, "presence_unavailable", "Presence index is not published")
    resolved = db.resolve()
    state: _Reader | None = getattr(_LOCAL, "reader", None)
    if state is not None and state.path == resolved and state.mtime == resolved.stat().st_mtime:
        return state
    close_presence()
    _LOCAL.reader = _open(resolved)
    return _LOCAL.reader


def scope_mask(conn: sqlite3.Connection, all_mask: int,
               include: Sequence[int], exclude: Sequence[int]) -> int:
    """Bitset of annotations in the taxon scope (lineage OR for include, minus exclude)."""
    def union(taxids: Sequence[int]) -> int:
        mask = 0
        if not taxids:
            return mask
        ph = ",".join("?" for _ in taxids)
        for (o,) in conn.execute(f"SELECT DISTINCT ord FROM lineage WHERE taxid IN ({ph})",
                                 list(taxids)):
            mask |= 1 << int(o)
        return mask

    mask = union(include) if include else all_mask
    if exclude:
        mask &= ~union(exclude)
    return mask & all_mask
