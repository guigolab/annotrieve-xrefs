"""
Manifest: the record of which annotations are indexed and under which ordinal.

``manifest.sqlite``. Ordinals are assigned
once and never reused. This module does **not** reconcile against Annotrieve
(decision D6, postponed): it registers the shards found on disk and stores the
catalogue metadata loaded by ``gene_corpus.sql.catalogue``.

Work is split so it can run in parallel without SQLite write contention:

* ``work_shard`` (any process): read the shard meta, classify the source, write the
  sidecar when needed. Touches no manifest. Returns a plain dict.
* ``record_result`` (the single manifest-owning process): register the annotation
  and store the outcome.
"""
from __future__ import annotations

import sqlite3
from collections.abc import Mapping
from pathlib import Path

from gene_corpus.sql.per_gff import SHARD_META_KEYS
from gene_corpus.sql.sidecar import (
    SHARD_NAME,
    SidecarError,
    is_stale,
    sidecar_path_for,
    write_sidecar,
)

MANIFEST_VERSION = "2"

STATUS_INDEXED = "indexed"
STATUS_PENDING = "pending"
STATUS_FAILED = "failed"
STATUS_SKIPPED_COMMUNITY = "skipped_community"
STATUS_REMOVED = "removed"

# Allow-list: only these ``source_database`` values (casefolded) are ever indexed.
# Community annotations (community_annotations.tsv of the annotation tracker, already
# skipped upstream by the import) are recorded as skipped if they ever show up: any
# value containing "community" (for example "Community" or "CommunityRegistry").
# Every other unknown value is recorded as ``failed`` so it is visible, never indexed.
INDEXED_SOURCES = frozenset({"ensembl", "genbank", "refseq"})
COMMUNITY_MARKER = "community"

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
    n_genes_named INTEGER,
    n_name_rows INTEGER,
    in_catalogue INTEGER NOT NULL DEFAULT 0,
    status TEXT NOT NULL,
    error TEXT
);
CREATE INDEX IF NOT EXISTS annotation_status ON annotation(status);

-- One row per annotation in the live Annotrieve catalogue (report endpoint, D8).
CREATE TABLE IF NOT EXISTS catalogue (
    annotation_id TEXT PRIMARY KEY,
    assembly_accession TEXT,
    assembly_name TEXT,
    organism_name TEXT,
    taxid INTEGER,
    database TEXT,
    provider TEXT,
    source_url TEXT,
    bgzip_path TEXT,
    release_date TEXT,
    taxon_lineage TEXT              -- ';'-joined taxids as served by the report
) WITHOUT ROWID;

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
        "INSERT INTO meta (key, value) VALUES ('manifest_version', ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (MANIFEST_VERSION,),
    )
    conn.commit()
    return conn


def classify_source(source_database: str | None) -> str | None:
    """Return ``STATUS_SKIPPED_COMMUNITY`` / ``STATUS_FAILED`` for unindexed sources, else None."""
    key = (source_database or "").strip().casefold()
    if key in INDEXED_SOURCES:
        return None
    if COMMUNITY_MARKER in key:
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


# --------------------------------------------------------------------------
# worker side (no manifest access)
# --------------------------------------------------------------------------
def work_shard(per_gff_root: str | Path, annotation_dir: str, force: bool = False) -> dict:
    """
    Inspect one shard and write its sidecar when it is allowed and stale.

    Returns a picklable dict: ``dir``, ``status`` (a STATUS_* value), ``meta`` (shard
    meta, when readable), ``error`` and, when a sidecar was written, ``info`` (counts and
    fingerprints). Never raises for a bad shard.
    """
    shard_path = Path(per_gff_root) / annotation_dir / SHARD_NAME
    result: dict = {"dir": annotation_dir, "status": STATUS_FAILED, "meta": None,
                    "error": None, "info": None, "written": False}
    try:
        meta = _read_shard_meta_ro(shard_path)
        missing = [k for k in SHARD_META_KEYS if k not in meta]
        if missing or not meta.get("annotation_id"):
            raise SidecarError(f"shard meta incomplete (missing {missing})")
    except Exception as exc:  # unreadable shard: nothing to register under
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result
    result["meta"] = meta

    skip = classify_source(meta.get("source_database"))
    if skip == STATUS_SKIPPED_COMMUNITY:
        result["status"] = skip
        return result
    if skip == STATUS_FAILED:
        result["error"] = f"unknown source_database {meta.get('source_database')!r}"
        return result

    try:
        sidecar = sidecar_path_for(shard_path)
        if force or is_stale(shard_path, sidecar):
            info = write_sidecar(shard_path, sidecar)
            result["written"] = True
            result["info"] = {
                "shard_fingerprint": info.shard_fingerprint,
                "sidecar_fingerprint": info.sidecar_fingerprint,
                "n_genes": info.n_genes,
                "n_genes_named": info.n_genes_named,
                "n_name_rows": info.n_name_rows,
                "size_bytes": info.size_bytes,
            }
        result["status"] = STATUS_INDEXED
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


# --------------------------------------------------------------------------
# manifest side
# --------------------------------------------------------------------------
def record_result(conn: sqlite3.Connection, result: Mapping) -> str:
    """Store one ``work_shard`` outcome. Does not commit; returns the final status."""
    meta = result.get("meta")
    if not meta:
        annotation_dir = result["dir"]
        conn.execute(
            "INSERT INTO annotation (ord, annotation_id, source_database, status, error) "
            "VALUES ((SELECT COALESCE(MAX(ord), 0) + 1 FROM annotation), ?, '', ?, ?) "
            "ON CONFLICT(annotation_id) DO UPDATE SET status=excluded.status, error=excluded.error",
            (annotation_dir, STATUS_FAILED, result.get("error")),
        )
        return STATUS_FAILED
    ord_ = register_annotation(conn, meta)
    info = result.get("info")
    if info:
        conn.execute(
            "UPDATE annotation SET shard_fingerprint=?, sidecar_fingerprint=?, n_genes=?, "
            "n_genes_named=?, n_name_rows=? WHERE ord=?",
            (info["shard_fingerprint"], info["sidecar_fingerprint"], info["n_genes"],
             info["n_genes_named"], info["n_name_rows"], ord_),
        )
    conn.execute(
        "UPDATE annotation SET status=?, error=? WHERE ord=?",
        (result["status"], result.get("error"), ord_),
    )
    return result["status"]
