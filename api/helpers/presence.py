"""Read-only thread-local connection to ``presence.sqlite`` and taxon-scope bitsets."""
from __future__ import annotations

import fcntl
import logging
import os
import shutil
import sqlite3
import threading
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from helpers.errors import http_error
from helpers.postings import iter_bits  # noqa: F401  (re-export for services)
from settings import presence_dir, presence_ram_dir, settings

PRESENCE_FILE = "presence.sqlite"
_LOCAL = threading.local()
_SOFT_HEAP_LIMIT = 64 * 1024 * 1024
_LOG = logging.getLogger(__name__)
_RAM_FAILED: set[Path] = set()  # versions that could not be copied to RAM: do not retry per request
_RAM_LOCK = threading.Lock()


@dataclass
class _Reader:
    path: Path
    mtime: float
    conn: sqlite3.Connection
    all_mask: int
    data_version: str
    ann_info: dict | None = None
    ann_paths: dict | None = None
    db_path: Path | None = None  # file actually opened (RAM copy or `path`)
    checked: float = 0.0  # monotonic time of the last check of the published link


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


def _same_file(a: os.stat_result, b: os.stat_result) -> bool:
    return a.st_size == b.st_size and a.st_mtime_ns == b.st_mtime_ns


def _prune_ram_dir(ram: Path, keep: str) -> None:
    """Keep the version in use and the newest other one (a worker may still be switching)."""
    others = sorted((d for d in ram.iterdir() if d.is_dir() and d.name != keep),
                    key=lambda d: d.stat().st_mtime, reverse=True)
    for d in others[1:]:
        shutil.rmtree(d, ignore_errors=True)


def ram_copy_path(resolved: Path) -> Path | None:
    """
    Path of the tmpfs copy of ``resolved`` (a ``presence.sqlite``), creating it when needed.

    One copy is shared by every worker (a file lock makes the first one copy, the rest wait).
    Returns None, and the caller reads ``resolved`` in place, when RAM mode is off, the file is
    larger than PRESENCE_RAM_MAX_MB, the tmpfs is too small, or the copy fails.
    """
    ram = presence_ram_dir()
    if ram is None or resolved in _RAM_FAILED:
        return None
    target = ram / resolved.parent.name / PRESENCE_FILE
    try:
        src = resolved.stat()
        if target.is_file() and _same_file(src, target.stat()):
            return target
        if src.st_size > settings.PRESENCE_RAM_MAX_MB * 1024 * 1024:
            raise OSError(f"{src.st_size} bytes is larger than PRESENCE_RAM_MAX_MB="
                          f"{settings.PRESENCE_RAM_MAX_MB}")
        ram.mkdir(parents=True, exist_ok=True)
        with open(ram / ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)  # released when the file closes
            if target.is_file() and _same_file(src, target.stat()):
                return target  # another worker copied it while we waited
            if shutil.disk_usage(ram).free < src.st_size + 16 * 1024 * 1024:
                raise OSError(f"not enough free space in {ram} for {src.st_size} bytes")
            target.parent.mkdir(parents=True, exist_ok=True)
            tmp = target.with_name(f"{PRESENCE_FILE}.{os.getpid()}.tmp")
            try:
                shutil.copyfile(resolved, tmp)
                os.utime(tmp, ns=(src.st_atime_ns, src.st_mtime_ns))
                if not _same_file(src, tmp.stat()):
                    raise OSError("copy has a different size than the source")
                os.replace(tmp, target)
            finally:
                tmp.unlink(missing_ok=True)
            _prune_ram_dir(ram, target.parent.name)
        _LOG.info("presence index %s copied to RAM (%s, %.0f MB)", resolved.parent.name, target,
                  src.st_size / 1e6)
        return target
    except OSError as exc:
        with _RAM_LOCK:
            _RAM_FAILED.add(resolved)
        _LOG.warning("presence index stays on disk (%s): %s", resolved, exc)
        return None


def _open(resolved: Path) -> _Reader:
    ram_file = ram_copy_path(resolved)
    db_file = ram_file or resolved
    conn = sqlite3.connect(f"file:{db_file}?mode=ro&immutable=1", uri=True,
                           check_same_thread=True)
    conn.row_factory = sqlite3.Row
    if ram_file is not None:
        # Pages come straight from the shared tmpfs file: no big private cache needed.
        conn.execute(f"PRAGMA cache_size=-{settings.PRESENCE_RAM_CACHE_KIB}")
        conn.execute(f"PRAGMA mmap_size={max(ram_file.stat().st_size, 1)}")
    else:
        conn.execute(f"PRAGMA cache_size={settings.READ_CACHE_SIZE}")
        conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute(f"PRAGMA soft_heap_limit={_SOFT_HEAP_LIMIT}")
    conn.execute("PRAGMA query_only=ON")
    mask = 0
    for (o,) in conn.execute("SELECT ord FROM annotation"):
        mask |= 1 << int(o)
    meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
    return _Reader(resolved, resolved.stat().st_mtime, conn, mask,
                   str(meta.get("data_version", "")), db_path=db_file,
                   checked=time.monotonic())


def presence_reader() -> _Reader:
    """This thread's reader; reopens when ``current`` points at a new version."""
    state: _Reader | None = getattr(_LOCAL, "reader", None)
    recheck = settings.PRESENCE_RECHECK_S
    if state is not None and recheck > 0 and time.monotonic() - state.checked < recheck:
        return state  # no network-drive calls between checks
    db = presence_db_path()
    if db is None or not db.is_file():
        raise http_error(503, "presence_unavailable", "Presence index is not published")
    resolved = db.resolve()
    if state is not None and state.path == resolved and state.mtime == resolved.stat().st_mtime:
        state.checked = time.monotonic()
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

    pres = connect(reader.db_path or reader.path)
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
