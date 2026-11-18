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
    ann_info: dict | None = None
    ann_paths: dict | None = None


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


LOCI_FILE = "loci.sqlite"


def loci_db_path() -> Path | None:
    root = presence_dir()
    return None if root is None else root / "current" / LOCI_FILE


def open_request_connections() -> tuple[sqlite3.Connection, sqlite3.Connection, str]:
    """
    Dedicated read-only connections ``(presence, loci, data_version)`` for one /loci request.

    Streaming responses are iterated from different worker threads, so these connections are
    not thread-bound; each belongs to one request and is used sequentially. The caller must
    close them. Raises 503 when loci.sqlite is missing or built for another presence version.
    """
    reader = presence_reader()
    path = loci_db_path()
    if path is None or not path.is_file():
        raise http_error(503, "loci_unavailable", "Loci index is not published")
    # Resolve both through the same version directory so a switch of `current` cannot split them.
    version_dir = reader.path.parent
    loci_path = version_dir / LOCI_FILE
    if not loci_path.is_file():
        raise http_error(503, "loci_unavailable", "Loci index is not published for this version")

    def connect(p: Path) -> sqlite3.Connection:
        conn = sqlite3.connect(f"file:{p}?mode=ro&immutable=1", uri=True, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        conn.execute(f"PRAGMA cache_size=-{settings.LOCI_CACHE_KIB}")
        conn.execute("PRAGMA temp_store=FILE")
        conn.execute("PRAGMA mmap_size=0")
        conn.execute("PRAGMA query_only=ON")
        return conn

    pres = connect(reader.path)
    loci = connect(loci_path)
    try:
        v = loci.execute("SELECT value FROM meta WHERE key='data_version'").fetchone()
        if v is None or v[0] != reader.data_version:
            raise http_error(503, "loci_version_mismatch",
                             "loci.sqlite was built for a different presence version")
    except Exception:
        pres.close()
        loci.close()
        raise
    return pres, loci, reader.data_version


def annotation_info(reader: "_Reader") -> dict[int, tuple]:
    """ord -> (annotation_id, organism_name, assembly_accession, source); loaded once per reader."""
    if reader.ann_info is None:
        reader.ann_info = {
            int(r["ord"]): (r["annotation_id"], r["organism_name"], r["assembly_accession"],
                            int(r["source"]))
            for r in reader.conn.execute(
                "SELECT ord, annotation_id, organism_name, assembly_accession, source FROM annotation")
        }
    return reader.ann_info


REF_FILE = "loci_ref.sqlite"


def annotation_paths(reader: "_Reader") -> dict[int, tuple[int, str]]:
    """ord -> (taxid, gff_path); loaded once per reader (used to find GFF files)."""
    if reader.ann_paths is None:
        reader.ann_paths = {
            int(r["ord"]): (int(r["taxid"]), r["gff_path"] or "")
            for r in reader.conn.execute("SELECT ord, taxid, gff_path FROM annotation")
        }
    return reader.ann_paths


def open_ref_connection() -> sqlite3.Connection:
    """Read-only connection to ``loci_ref.sqlite`` (``locus_ref`` -> locus lookup) for one request.

    Raises 503 when the file is missing (build it with ``serve_build loci-ref``) or was built for
    another presence version. The caller closes it.
    """
    reader = presence_reader()
    path = reader.path.parent / REF_FILE
    if not path.is_file():
        raise http_error(503, "locus_ref_index_unavailable",
                         "locus_ref lookup index is not published (run `serve_build loci-ref`)")
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA query_only=ON")
    row = conn.execute("SELECT value FROM meta WHERE key='data_version'").fetchone()
    if row is None or row[0] != reader.data_version:
        conn.close()
        raise http_error(503, "loci_version_mismatch",
                         "loci_ref.sqlite was built for a different presence version")
    return conn
