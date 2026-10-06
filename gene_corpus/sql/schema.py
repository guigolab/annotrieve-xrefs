"""
Global SQLite schema for annotrieve-xrefs serve (sidecar locked).

Tables: annotation, annotation_lineage, gene_hit, xref_meta, namespace_stats.
``xref_hit`` is not written.

Incremental attach (``sql.attach`` / ``gene_corpus.sync``):
  1. INSERT annotation row with next annotation_key
  2. INSERT gene_hit from shard gene_xref ⋈ gene (Tier A)
  3. Upsert xref_meta: n_annotations += 1, n_loci += gene_count
  4. Rebuild namespace_stats from xref_meta (or bump when accession is new)
  5. Insert annotation_lineage for the new key
  6. Readers reopen on mtime (immutable=1 while writers are offline)
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

SCHEMA_VERSION = 1

# Negative cache_size is KiB.
_DEFAULT_BUILD_CACHE_SIZE = -131072  # 128 MiB
# Cap SQLite heap so large sorts spill to temp files (sidecar RAM budget).
_SOFT_HEAP_LIMIT = 2 * 1024 * 1024 * 1024  # 2 GiB

_TABLES_SQL = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS annotation (
    annotation_key INTEGER PRIMARY KEY,
    annotation_id TEXT UNIQUE NOT NULL,
    taxid INTEGER NOT NULL,
    assembly_accession TEXT NOT NULL,
    organism_name TEXT,
    source_database TEXT NOT NULL,
    source_provider TEXT,
    profile_id TEXT NOT NULL,
    gff_path TEXT
);

CREATE TABLE IF NOT EXISTS annotation_lineage (
    annotation_key INTEGER NOT NULL REFERENCES annotation(annotation_key),
    taxid INTEGER NOT NULL,
    PRIMARY KEY (annotation_key, taxid)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS gene_hit (
    namespace TEXT NOT NULL,
    accession TEXT NOT NULL,
    annotation_key INTEGER NOT NULL REFERENCES annotation(annotation_key),
    local_id INTEGER NOT NULL,
    seqid TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand INTEGER NOT NULL,
    feature_type TEXT NOT NULL,
    biotype TEXT,
    primary_name TEXT,
    PRIMARY KEY (namespace, accession, annotation_key, local_id)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS xref_meta (
    namespace TEXT NOT NULL,
    accession TEXT NOT NULL,
    n_annotations INTEGER NOT NULL,
    n_loci INTEGER NOT NULL,
    PRIMARY KEY (namespace, accession)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS namespace_stats (
    namespace TEXT PRIMARY KEY,
    accession_count INTEGER NOT NULL
) WITHOUT ROWID;
"""

_INDEXES_SQL = """
CREATE INDEX IF NOT EXISTS annotation_taxid_idx
    ON annotation (taxid);

CREATE INDEX IF NOT EXISTS annotation_lineage_taxid
    ON annotation_lineage (taxid, annotation_key);

CREATE INDEX IF NOT EXISTS xref_meta_ns_loci
    ON xref_meta (namespace, n_loci DESC, accession);

CREATE INDEX IF NOT EXISTS xref_meta_ns_count
    ON xref_meta (namespace, n_annotations DESC, accession);
"""


def connect_for_build(path: str | Path) -> sqlite3.Connection:
    """Open SQLite for build (WAL, modest cache, soft heap cap)."""
    conn = sqlite3.connect(str(path))
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute(f"PRAGMA cache_size={_DEFAULT_BUILD_CACHE_SIZE}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute(f"PRAGMA soft_heap_limit={_SOFT_HEAP_LIMIT}")
    return conn


def init_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(_TABLES_SQL)
    conn.execute(
        "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
        ("schema_version", str(SCHEMA_VERSION)),
    )
    conn.commit()


def create_indexes(conn: sqlite3.Connection) -> None:
    conn.executescript(_INDEXES_SQL)
    conn.execute("ANALYZE")
    conn.commit()


def seed_annotations(
    conn: sqlite3.Connection,
    rows: list[dict],
) -> None:
    """
    Insert annotation dimension rows.

    Each dict: annotation_key, annotation_id, taxid, assembly_accession,
    organism_name, source_database, source_provider, profile_id, gff_path.
    """
    conn.executemany(
        """
        INSERT OR REPLACE INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id, gff_path
        ) VALUES (
            :annotation_key, :annotation_id, :taxid, :assembly_accession,
            :organism_name, :source_database, :source_provider, :profile_id, :gff_path
        )
        """,
        rows,
    )
    conn.commit()
