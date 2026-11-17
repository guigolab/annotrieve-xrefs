"""
Manifest: the record of which annotations are indexed and under which ordinal.

``manifest.sqlite`` (see ``PLAN-build-tables.md`` section 4). Ordinals are assigned
once and never reused. This module does **not** reconcile against Annotrieve
(decision D6, postponed): it only registers shards found on disk.
"""
from __future__ import annotations

import sqlite3
from collections.abc import Mapping
from pathlib import Path

from gene_corpus.sql.per_gff import SHARD_META_KEYS
from gene_corpus.sql.sidecar import (
    SHARD_NAME,
    SidecarError,
    SidecarInfo,
    is_stale,
    sidecar_path_for,
    write_sidecar,
)

MANIFEST_VERSION = "1"

STATUS_INDEXED = "indexed"
STATUS_PENDING = "pending"
STATUS_FAILED = "failed"
STATUS_SKIPPED_COMMUNITY = "skipped_community"
STATUS_REMOVED = "removed"

# source_database values that are indexed. Compared casefolded. Anything else is
# recorded with its own status and not indexed (D9: confirm the real spellings in M2).
INDEXED_SOURCES = frozenset({"ensembl", "genbank", "refseq"})
COMMUNITY_SOURCE = "community"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS annotation (
    ord INTEGER PRIMARY KEY,
    annotation_id TEXT NOT NULL UNIQUE,
    taxid INTEGER,
    organism_name TEXT,
    assembly_accession TEXT,
    assembly_name TEXT,
    source_database TEXT NOT NULL,
    source_provider TEXT,
    release_date TEXT,
    gff_path TEXT,
    shard_fingerprint TEXT,
    sidecar_fingerprint TEXT,
    n_genes INTEGER,
    n_name_rows INTEGER,
    status TEXT NOT NULL,
    error TEXT
);
CREATE INDEX IF NOT EXISTS annotation_status ON annotation(status);

CREATE TABLE IF NOT EXISTS lineage (
    taxid INTEGER NOT NULL,
    ord INTEGER NOT NULL,
    PRIMARY KEY (taxid, ord)
) WITHOUT ROWID;
"""


def open_manifest(path: str | Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(path))
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.executescript(_SCHEMA)
    conn.execute(
        "INSERT OR IGNORE INTO meta (key, value) VALUES ('manifest_version', ?)",
        (MANIFEST_VERSION,),
    )
    conn.commit()
    return conn


def classify_source(source_database: str | None) -> str | None:
    """Return ``STATUS_SKIPPED_COMMUNITY`` / ``STATUS_FAILED`` for unindexed sources, else None."""
    key = (source_database or "").strip().casefold()
    if key in INDEXED_SOURCES:
        return None
    if key == COMMUNITY_SOURCE:
        return STATUS_SKIPPED_COMMUNITY
    return STATUS_FAILED


def _int_or_none(value: str | None) -> int | None:
    try:
        return int(value) if value not in (None, "") else None
    except ValueError:
        return None


def register_annotation(conn: sqlite3.Connection, meta: Mapping[str, str]) -> int:
    """Insert or update the manifest row for a shard's meta; return its ``ord``."""
    annotation_id = meta["annotation_id"]
    row = conn.execute(
        "SELECT ord FROM annotation WHERE annotation_id = ?", (annotation_id,)
    ).fetchone()
    values = (
        _int_or_none(meta.get("taxid")),
        meta.get("organism_name"),
        meta.get("assembly_accession"),
        meta.get("source_database") or "",
        meta.get("source_provider"),
        meta.get("gff_path"),
    )
    if row:
        conn.execute(
            "UPDATE annotation SET taxid=?, organism_name=?, assembly_accession=?, "
            "source_database=?, source_provider=?, gff_path=? WHERE ord=?",
            (*values, row[0]),
        )
        return int(row[0])
    nxt = conn.execute("SELECT COALESCE(MAX(ord), 0) + 1 FROM annotation").fetchone()[0]
    conn.execute(
        "INSERT INTO annotation (ord, annotation_id, taxid, organism_name, "
        "assembly_accession, source_database, source_provider, gff_path, status) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (nxt, annotation_id, *values, STATUS_PENDING),
    )
    return int(nxt)


def _read_shard_meta_ro(shard_path: Path) -> dict[str, str]:
    conn = sqlite3.connect(f"file:{shard_path}?mode=ro", uri=True)
    try:
        return {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
    finally:
        conn.close()


def _set_status(conn, ord_: int, status: str, error: str | None = None) -> None:
    conn.execute("UPDATE annotation SET status=?, error=? WHERE ord=?", (status, error, ord_))


def index_shard(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    annotation_dir: str,
    *,
    force: bool = False,
) -> str:
    """
    Register one shard and (re)write its sidecar when needed. Returns the final status.

    Never raises for a bad shard: the failure is recorded as ``failed`` with a message.
    Commits after each shard so a crash loses at most the shard in flight.
    """
    shard_path = Path(per_gff_root) / annotation_dir / SHARD_NAME
    try:
        meta = _read_shard_meta_ro(shard_path)
        missing = [k for k in SHARD_META_KEYS if k not in meta]
        if missing or not meta.get("annotation_id"):
            raise SidecarError(f"shard meta incomplete (missing {missing})")
    except Exception as exc:  # unreadable shard: nothing to register under
        return _record_unreadable(conn, annotation_dir, exc)

    ord_ = register_annotation(conn, meta)
    skip = classify_source(meta.get("source_database"))
    if skip == STATUS_SKIPPED_COMMUNITY:
        _set_status(conn, ord_, skip)
        conn.commit()
        return skip
    if skip == STATUS_FAILED:
        _set_status(conn, ord_, STATUS_FAILED, f"unknown source_database {meta.get('source_database')!r}")
        conn.commit()
        return STATUS_FAILED

    try:
        sidecar = sidecar_path_for(shard_path)
        if force or is_stale(shard_path, sidecar):
            info: SidecarInfo = write_sidecar(shard_path, sidecar)
            conn.execute(
                "UPDATE annotation SET shard_fingerprint=?, sidecar_fingerprint=?, "
                "n_genes=?, n_name_rows=? WHERE ord=?",
                (info.shard_fingerprint, info.sidecar_fingerprint, info.n_genes,
                 info.n_name_rows, ord_),
            )
        _set_status(conn, ord_, STATUS_INDEXED)
    except Exception as exc:
        _set_status(conn, ord_, STATUS_FAILED, f"{type(exc).__name__}: {exc}")
    conn.commit()
    return conn.execute("SELECT status FROM annotation WHERE ord=?", (ord_,)).fetchone()[0]


def _record_unreadable(conn: sqlite3.Connection, annotation_dir: str, exc: Exception) -> str:
    """A shard we cannot even read meta from is recorded by directory name."""
    conn.execute(
        "INSERT INTO annotation (ord, annotation_id, source_database, status, error) "
        "VALUES ((SELECT COALESCE(MAX(ord), 0) + 1 FROM annotation), ?, '', ?, ?) "
        "ON CONFLICT(annotation_id) DO UPDATE SET status=excluded.status, error=excluded.error",
        (annotation_dir, STATUS_FAILED, f"{type(exc).__name__}: {exc}"),
    )
    conn.commit()
    return STATUS_FAILED
