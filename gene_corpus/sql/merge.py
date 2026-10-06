"""Merge per-GFF shards into global gene_hit + xref_meta (+ lineage)."""
from __future__ import annotations

import sqlite3
import sys
from collections.abc import Mapping
from pathlib import Path

from gene_corpus.namespaces.tiers import TIER_A
from gene_corpus.sql.lineage import fill_annotation_lineage
from gene_corpus.sql.per_gff import shard_sqlite_path
from gene_corpus.sql.schema import create_indexes


def _tier_a_placeholders() -> tuple[str, list[str]]:
    ns_list = sorted(TIER_A)
    return ", ".join("?" for _ in ns_list), ns_list


def merge_shards_into_db(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    *,
    taxonomy_parents: Mapping[int, int | None] | None = None,
) -> tuple[int, int, int]:
    """
    Stream each shard into ``gene_hit`` and upsert ``xref_meta``.

    One shard attached at a time. Does not ``GROUP BY gene_hit``.
    When ``taxonomy_parents`` is set, fills ``annotation_lineage`` after
    ``namespace_stats``.

    Raises ``FileNotFoundError`` if any seeded annotation lacks a shard.
    Attach/stream failures abort the merge (no orphan annotation rows with
    hits partially applied for that shard).

    Returns ``(n_gene_hit, n_xref_meta, n_lineage)``.
    """
    conn.execute("DELETE FROM gene_hit")
    conn.execute("DELETE FROM xref_meta")
    conn.execute("DELETE FROM namespace_stats")
    conn.execute("DELETE FROM annotation_lineage")
    conn.commit()

    # Materialize first: ATTACH/INSERT cannot run while a SELECT cursor is
    # open on the same connection. Annotation dimension is small (key + id).
    ann_rows = conn.execute(
        "SELECT annotation_key, annotation_id FROM annotation "
        "ORDER BY annotation_key"
    ).fetchall()

    total = len(ann_rows)
    progress_every = 50
    done = 0
    ns_sql, ns_list = _tier_a_placeholders()
    print(
        f"merge attach: starting for {total} annotation(s)",
        file=sys.stderr,
        flush=True,
    )

    for annotation_key, annotation_id in ann_rows:
        shard = shard_sqlite_path(per_gff_root, str(annotation_id))
        if not shard.is_file():
            raise FileNotFoundError(
                f"missing shard for {annotation_id}: {shard}"
            )

        # Absolute path: ATTACH is relative to the process cwd, not the DB.
        try:
            conn.execute(
                "ATTACH DATABASE ? AS shard", (str(shard.resolve()),)
            )
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
                    [int(annotation_key), *ns_list],
                )
                # WHERE true: SQLite parse quirk for INSERT…SELECT…ON CONFLICT
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
                # Commit before DETACH — SQLite holds an attach lock until commit.
                conn.commit()
            except Exception:
                conn.rollback()
                raise
        finally:
            conn.execute("DETACH DATABASE shard")

        done += 1
        if done % progress_every == 0 or done == total:
            print(
                f"merge progress: {done}/{total}",
                file=sys.stderr,
                flush=True,
            )

    print(
        "merge: building namespace_stats…",
        file=sys.stderr,
        flush=True,
    )
    conn.execute(
        """
        INSERT INTO namespace_stats (namespace, accession_count)
        SELECT namespace, COUNT(*)
        FROM xref_meta
        GROUP BY namespace
        """
    )
    conn.commit()

    n_lineage = 0
    if taxonomy_parents is not None:
        print("merge: filling annotation_lineage…", file=sys.stderr, flush=True)
        n_lineage = fill_annotation_lineage(conn, taxonomy_parents)

    n_hit = int(conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0])
    n_meta = int(conn.execute("SELECT COUNT(*) FROM xref_meta").fetchone()[0])

    print("merge: ANALYZE / indexes…", file=sys.stderr, flush=True)
    create_indexes(conn)
    return n_hit, n_meta, n_lineage
