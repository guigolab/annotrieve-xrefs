"""Merge per-GFF shards into global gene_hit + xref_meta (+ lineage)."""
from __future__ import annotations

import sqlite3
import sys
from collections.abc import Mapping
from pathlib import Path

from gene_corpus.namespaces.tiers import TIER_A
from gene_corpus.sql.lineage import fill_annotation_lineage
from gene_corpus.sql.per_gff import shard_sqlite_path
from gene_corpus.sql.schema import checkpoint_wal, create_indexes


def _tier_a_placeholders() -> tuple[str, list[str]]:
    ns_list = sorted(TIER_A)
    return ", ".join("?" for _ in ns_list), ns_list


def _is_gene_hit_duplicate(exc: sqlite3.IntegrityError) -> bool:
    """True only for a UNIQUE/PK conflict on ``gene_hit`` (key already merged)."""
    return str(exc).startswith("UNIQUE constraint failed: gene_hit.")


def _merge_done_count(conn: sqlite3.Connection) -> int:
    row = conn.execute("SELECT COUNT(*) FROM merge_done").fetchone()
    return int(row[0]) if row else 0


def _mark_merge_done_prefix(conn: sqlite3.Connection, through_key: int) -> int:
    """INSERT OR IGNORE keys 1..through_key into merge_done. Returns rows attempted."""
    if through_key < 1:
        return 0
    conn.executemany(
        "INSERT OR IGNORE INTO merge_done(annotation_key) VALUES (?)",
        [(k,) for k in range(1, through_key + 1)],
    )
    conn.commit()
    return through_key


def _bootstrap_merge_done(
    conn: sqlite3.Connection,
    *,
    resume_after_key: int | None,
) -> None:
    """
    Seed merge_done for a legacy interrupted DB (no progress rows yet).

    Prefer resume_after_key (O(1)). Otherwise MAX(gene_hit.annotation_key)
    (correct contiguous prefix; may be slow on large NFS DBs).
    """
    if resume_after_key is not None:
        if resume_after_key < 0:
            raise ValueError("resume_after_key must be >= 0")
        if resume_after_key >= 1:
            n = _mark_merge_done_prefix(conn, resume_after_key)
            print(
                f"merge resume: seeded merge_done for keys 1..{n} "
                f"(--resume-after-key)",
                file=sys.stderr,
                flush=True,
            )
        return

    if _merge_done_count(conn) > 0:
        return

    print(
        "merge resume: merge_done empty; scanning MAX(gene_hit.annotation_key) "
        "(slow on large NFS DBs; prefer --resume-after-key)",
        file=sys.stderr,
        flush=True,
    )
    row = conn.execute("SELECT MAX(annotation_key) FROM gene_hit").fetchone()
    max_key = row[0] if row else None
    if max_key is None:
        print(
            "merge resume: no gene_hit rows; starting from key 1",
            file=sys.stderr,
            flush=True,
        )
        return
    k = int(max_key)
    _mark_merge_done_prefix(conn, k)
    print(
        f"merge resume: seeded merge_done for keys 1..{k} from gene_hit MAX",
        file=sys.stderr,
        flush=True,
    )


def merge_shards_into_db(
    conn: sqlite3.Connection,
    per_gff_root: Path,
    *,
    taxonomy_parents: Mapping[int, int | None] | None = None,
    resume: bool = False,
    resume_after_key: int | None = None,
    checkpoint_every: int = 50,
) -> tuple[int, int, int]:
    """
    Stream each shard into ``gene_hit`` and upsert ``xref_meta``.

    One shard attached at a time. Does not ``GROUP BY gene_hit``.
    When ``taxonomy_parents`` is set, fills ``annotation_lineage`` after
    ``namespace_stats``.

    When ``resume`` is True, keeps existing ``gene_hit`` / ``xref_meta`` /
    ``merge_done`` and skips keys already in ``merge_done``. Optional
    ``resume_after_key`` seeds ``merge_done`` for ``1..N`` (legacy cutoff).

    Raises ``FileNotFoundError`` if any seeded annotation lacks a shard.
    Attach/stream failures abort the merge (no orphan annotation rows with
    hits partially applied for that shard), except UNIQUE conflicts on
    resume which mark the key done and continue.

    Returns ``(n_gene_hit, n_xref_meta, n_lineage)``.
    """
    if resume_after_key is not None and not resume:
        raise ValueError("resume_after_key requires resume=True")

    if resume:
        conn.execute("DELETE FROM namespace_stats")
        conn.execute("DELETE FROM annotation_lineage")
        conn.commit()
        if resume_after_key is not None:
            _bootstrap_merge_done(conn, resume_after_key=resume_after_key)
        elif _merge_done_count(conn) == 0:
            _bootstrap_merge_done(conn, resume_after_key=None)
    else:
        conn.execute("DELETE FROM gene_hit")
        conn.execute("DELETE FROM xref_meta")
        conn.execute("DELETE FROM namespace_stats")
        conn.execute("DELETE FROM annotation_lineage")
        conn.execute("DELETE FROM merge_done")
        conn.commit()

    # Materialize first: ATTACH/INSERT cannot run while a SELECT cursor is
    # open on the same connection. Annotation dimension is small (key + id).
    ann_rows = conn.execute(
        "SELECT annotation_key, annotation_id FROM annotation "
        "ORDER BY annotation_key"
    ).fetchall()

    total = len(ann_rows)
    progress_every = 50
    if checkpoint_every < 1:
        checkpoint_every = progress_every
    done_keys = {
        int(r[0])
        for r in conn.execute("SELECT annotation_key FROM merge_done")
    }
    done = len(done_keys)
    ns_sql, ns_list = _tier_a_placeholders()
    print(
        f"merge attach: starting for {total} annotation(s) "
        f"(already_done={done} resume={resume})",
        file=sys.stderr,
        flush=True,
    )
    if done and (done % progress_every == 0 or done == total):
        print(
            f"merge progress: {done}/{total}",
            file=sys.stderr,
            flush=True,
        )

    newly_committed = 0
    for annotation_key, annotation_id in ann_rows:
        key_i = int(annotation_key)
        if key_i in done_keys:
            continue

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

        already_present = False
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
                    [key_i, *ns_list],
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
                conn.execute(
                    "INSERT OR IGNORE INTO merge_done(annotation_key) VALUES (?)",
                    (key_i,),
                )
                # Commit before DETACH — SQLite holds an attach lock until commit.
                conn.commit()
            except sqlite3.IntegrityError as exc:
                # Resume only: low --resume-after-key re-probes keys already
                # in gene_hit. Anything other than a gene_hit PK duplicate
                # (or any conflict on a fresh merge) must abort.
                conn.rollback()
                if not resume or not _is_gene_hit_duplicate(exc):
                    raise
                already_present = True
            except Exception:
                conn.rollback()
                raise
        finally:
            conn.execute("DETACH DATABASE shard")

        if already_present:
            print(
                f"merge resume: key {key_i} ({annotation_id}) already in "
                f"gene_hit; marking merge_done",
                file=sys.stderr,
                flush=True,
            )
            conn.execute(
                "INSERT OR IGNORE INTO merge_done(annotation_key) VALUES (?)",
                (key_i,),
            )
            conn.commit()

        done_keys.add(key_i)
        done += 1
        newly_committed += 1
        if done % progress_every == 0 or done == total:
            print(
                f"merge progress: {done}/{total}",
                file=sys.stderr,
                flush=True,
            )
        if newly_committed % checkpoint_every == 0:
            checkpoint_wal(conn)

    if newly_committed > 0 and newly_committed % checkpoint_every != 0:
        checkpoint_wal(conn)

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
