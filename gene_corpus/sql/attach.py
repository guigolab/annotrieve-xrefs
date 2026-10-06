"""Incremental attach of one per-GFF shard into the global gene_corpus DB."""
from __future__ import annotations

import sqlite3
import sys
from collections.abc import Mapping
from pathlib import Path

from gene_corpus.namespaces.tiers import TIER_A
from gene_corpus.sql.lineage import fill_lineage_for_annotation
from gene_corpus.sql.per_gff import (
    meta_complete,
    read_shard_meta,
    shard_sqlite_path,
)


def _tier_a_placeholders() -> tuple[str, list[str]]:
    ns_list = sorted(TIER_A)
    return ", ".join("?" for _ in ns_list), ns_list


def next_annotation_key(conn: sqlite3.Connection) -> int:
    """Return ``MAX(annotation_key) + 1``, or 1 when the table is empty."""
    row = conn.execute(
        "SELECT COALESCE(MAX(annotation_key), 0) FROM annotation"
    ).fetchone()
    return int(row[0]) + 1


def annotation_id_exists(conn: sqlite3.Connection, annotation_id: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM annotation WHERE annotation_id = ? LIMIT 1",
        (annotation_id,),
    ).fetchone()
    return row is not None


def existing_annotation_ids(conn: sqlite3.Connection) -> set[str]:
    return {
        str(r[0])
        for r in conn.execute("SELECT annotation_id FROM annotation")
    }


def seed_row_from_shard_meta(
    annotation_key: int,
    meta: Mapping[str, str],
) -> dict:
    """Build an annotation seed dict from shard meta + assigned key."""
    if not meta_complete(meta):
        raise ValueError("incomplete shard meta")
    ann_id = meta.get("annotation_id") or ""
    if not ann_id:
        raise ValueError("empty annotation_id in shard meta")
    return {
        "annotation_key": int(annotation_key),
        "annotation_id": ann_id,
        "taxid": int(meta["taxid"]),
        "assembly_accession": meta["assembly_accession"],
        "organism_name": meta["organism_name"] or None,
        "source_database": meta["source_database"],
        "source_provider": meta["source_provider"] or None,
        "profile_id": meta["profile_id"],
        "gff_path": meta["gff_path"] or None,
    }


def read_seed_row_from_shard(
    per_gff_root: Path,
    annotation_id: str,
    annotation_key: int,
) -> dict:
    """Open shard read-only and build a seed row for ``annotation_key``."""
    shard = shard_sqlite_path(per_gff_root, annotation_id)
    if not shard.is_file():
        raise FileNotFoundError(f"missing shard: {shard}")
    conn = sqlite3.connect(f"file:{shard.resolve()}?mode=ro", uri=True)
    try:
        meta = read_shard_meta(conn)
    finally:
        conn.close()
    seed = seed_row_from_shard_meta(annotation_key, meta)
    meta_id = str(seed["annotation_id"])
    if meta_id != annotation_id:
        raise ValueError(
            f"shard meta annotation_id={meta_id!r} does not match "
            f"directory id={annotation_id!r} ({shard})"
        )
    return seed


def rebuild_namespace_stats(conn: sqlite3.Connection) -> int:
    """Rebuild ``namespace_stats`` from ``xref_meta``. Returns row count."""
    conn.execute("DELETE FROM namespace_stats")
    conn.execute(
        """
        INSERT INTO namespace_stats (namespace, accession_count)
        SELECT namespace, COUNT(*)
        FROM xref_meta
        GROUP BY namespace
        """
    )
    conn.commit()
    return int(conn.execute("SELECT COUNT(*) FROM namespace_stats").fetchone()[0])


def attach_shard(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    seed_row: Mapping,
    *,
    taxonomy_parents: Mapping[int, int | None] | None = None,
) -> bool:
    """
    Insert one annotation and stream its shard into global tables.

    ``seed_row`` must include ``annotation_key`` and dimension fields.
    Idempotent: if ``annotation_id`` already exists, logs and returns False.

    When ``taxonomy_parents`` is set, fills lineage for this key only.
    Annotation / gene_hit / xref_meta / lineage commit as one transaction.
    Returns True when the shard was attached.
    """
    annotation_id = str(seed_row["annotation_id"])
    annotation_key = int(seed_row["annotation_key"])
    if annotation_id_exists(conn, annotation_id):
        print(
            f"attach skip: annotation_id already present ({annotation_id})",
            file=sys.stderr,
        )
        return False

    shard = shard_sqlite_path(per_gff_root, annotation_id)
    if not shard.is_file():
        raise FileNotFoundError(f"missing shard for {annotation_id}: {shard}")

    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id, gff_path
        ) VALUES (
            :annotation_key, :annotation_id, :taxid, :assembly_accession,
            :organism_name, :source_database, :source_provider, :profile_id,
            :gff_path
        )
        """,
        dict(seed_row),
    )

    ns_sql, ns_list = _tier_a_placeholders()
    try:
        conn.execute("ATTACH DATABASE ? AS shard", (str(shard.resolve()),))
    except Exception:
        conn.rollback()
        raise

    try:
        try:
            conn.execute(
                f"""
                INSERT INTO gene_hit (
                    namespace, accession, annotation_key, local_id,
                    seqid, start, end, strand,
                    feature_type, biotype, primary_name
                )
                SELECT
                    x.namespace, x.accession, ?, x.local_id,
                    g.seqid, g.start, g.end, g.strand,
                    g.feature_type, g.biotype, g.primary_name
                FROM shard.gene_xref x
                JOIN shard.gene g ON g.local_id = x.local_id
                WHERE x.namespace IN ({ns_sql})
                """,
                [annotation_key, *ns_list],
            )
            conn.execute(
                """
                INSERT INTO xref_meta (
                    namespace, accession, n_annotations, n_loci
                )
                SELECT namespace, accession, 1, gene_count
                FROM shard.tier_a_counts
                WHERE true
                ON CONFLICT(namespace, accession) DO UPDATE SET
                    n_annotations = n_annotations + 1,
                    n_loci = n_loci + excluded.n_loci
                """
            )
            if taxonomy_parents is not None:
                fill_lineage_for_annotation(
                    conn,
                    annotation_key,
                    int(seed_row["taxid"]),
                    taxonomy_parents,
                    commit=False,
                )
            # Commit before DETACH — SQLite holds an attach lock until commit.
            conn.commit()
        except Exception:
            conn.rollback()
            raise
    finally:
        conn.execute("DETACH DATABASE shard")

    return True
