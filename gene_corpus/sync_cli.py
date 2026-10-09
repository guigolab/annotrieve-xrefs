"""
Sync CLI: harvest new annotations and incrementally attach into gene_corpus.sqlite.

    python -m gene_corpus.sync \\
        --work-dir /data/gene_corpus \\
        --files-root /data/annotrieve/files \\
        --taxonomy-tsv /data/annotrieve/files/taxonomy/flattened-tree.tsv

TODO: Quote-wrapped symbol/alias accessions already present in
gene_corpus.sqlite are cleaned with ``python -m gene_corpus.repair`` for the
time being. Sync/attach does not rewrite historical keys; harvest
``normalize_symbol`` only prevents new wrapped accessions.
"""
from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from gene_corpus.catalog import DEFAULT_API_BASE, load_report_rows
from gene_corpus.harvest_cli import prepare_job, process_one
from gene_corpus.resolve import resolve_jobs
from gene_corpus.sql.attach import (
    attach_shard,
    existing_annotation_ids,
    next_annotation_key,
    read_seed_row_from_shard,
    rebuild_namespace_stats,
)
from gene_corpus.sql.lineage import (
    load_parent_map,
    refresh_lineage_if_taxonomy_changed,
)
from gene_corpus.sql.per_gff import shard_sqlite_path
from gene_corpus.sql.schema import connect_for_build

DEFAULT_TAXONOMY_TSV = Path(
    "/data/annotrieve/files/taxonomy/flattened-tree.tsv"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Harvest new annotations from the production report and "
            "incrementally attach them into an existing gene_corpus.sqlite. "
            "Refreshes annotation_lineage when the taxonomy TSV content changes."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="corpus root (gene_corpus.sqlite + per_gff/)",
    )
    parser.add_argument(
        "--files-root",
        type=Path,
        required=True,
        help="host mount of annotations volume (GFF paths)",
    )
    parser.add_argument(
        "--api-base",
        default=DEFAULT_API_BASE,
        help=f"Annotrieve API root (default: {DEFAULT_API_BASE})",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="local annotation report TSV (skip download)",
    )
    parser.add_argument(
        "--taxonomy-tsv",
        type=Path,
        default=DEFAULT_TAXONOMY_TSV,
        help=(
            "Annotrieve flattened taxonomy TSV on the host "
            f"(default: {DEFAULT_TAXONOMY_TSV})"
        ),
    )
    parser.add_argument(
        "--download-timeout",
        type=float,
        default=600.0,
        help="timeout seconds for annotation report download",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="harvest process pool size (1 or 2; default: 1)",
    )
    return parser


def _touch_db_for_readers(db_path: Path, conn: sqlite3.Connection) -> None:
    """Checkpoint WAL and bump mtime so API readers reopen."""
    try:
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    except sqlite3.Error as exc:
        print(f"sync warn: wal_checkpoint failed: {exc}", file=sys.stderr)
    now = time.time()
    os.utime(db_path, (now, now))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    per_gff_root = work_dir / "per_gff"
    db_path = work_dir / "gene_corpus.sqlite"
    taxonomy_tsv: Path = args.taxonomy_tsv
    files_root: Path = args.files_root

    if not db_path.is_file():
        print(
            f"gene_corpus.sqlite missing: {db_path} "
            "(run python -m gene_corpus.merge first)",
            file=sys.stderr,
        )
        return 1
    if not taxonomy_tsv.is_file():
        print(f"taxonomy TSV missing: {taxonomy_tsv}", file=sys.stderr)
        return 1
    if args.workers < 1 or args.workers > 2:
        print("--workers must be 1 or 2", file=sys.stderr)
        return 2

    per_gff_root.mkdir(parents=True, exist_ok=True)

    try:
        parents = load_parent_map(taxonomy_tsv)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(
        f"taxonomy: loaded {len(parents)} taxid→parent edges from {taxonomy_tsv}",
        file=sys.stderr,
        flush=True,
    )

    conn = connect_for_build(db_path)
    known_ids = existing_annotation_ids(conn)
    print(
        f"sync: {len(known_ids)} annotation(s) already in global DB",
        file=sys.stderr,
        flush=True,
    )

    try:
        rows = load_report_rows(
            api_base=None if args.report else args.api_base,
            report=args.report,
            annotation_ids=None,
            timeout=args.download_timeout,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        conn.close()
        return 1

    new_rows = [
        r for r in rows
        if r.annotation_id and r.annotation_id not in known_ids
    ]
    print(
        f"report rows: {len(rows)}; new (not in DB): {len(new_rows)}",
        file=sys.stderr,
        flush=True,
    )

    harvest_errors = 0
    harvested = 0
    skipped_existing = 0
    if new_rows:
        jobs, stats = resolve_jobs(
            new_rows,
            gff_mode="local",
            files_root=files_root,
        )
        print(
            f"resolve: ready={stats.ready} "
            f"skipped_missing={stats.skipped_missing} "
            f"skipped_no_url={stats.skipped_no_url} "
            f"skipped_no_profile={stats.skipped_no_profile}",
            file=sys.stderr,
        )
        prepared = [prepare_job(job) for job in jobs]
        total = len(prepared)
        if prepared:
            print(
                f"harvest: {total} job(s) for new annotation(s)",
                file=sys.stderr,
                flush=True,
            )

            def _on_complete(done: int) -> None:
                if done % 50 == 0 or done == total:
                    print(
                        f"harvest progress: {done}/{total} "
                        f"(ok={harvested} skipped={skipped_existing} "
                        f"errors={harvest_errors})",
                        file=sys.stderr,
                        flush=True,
                    )

            done = 0
            if args.workers == 1:
                for job in prepared:
                    try:
                        n_pairs, ann_id, skipped = process_one(
                            job, str(per_gff_root), force=False
                        )
                        if skipped:
                            skipped_existing += 1
                        else:
                            harvested += 1
                        if n_pairs == 0:
                            print(
                                f"harvest empty {ann_id}: 0 xref pairs",
                                file=sys.stderr,
                            )
                    except Exception as exc:  # noqa: BLE001
                        harvest_errors += 1
                        print(
                            f"harvest fail {job.annotation_id}: {exc}",
                            file=sys.stderr,
                        )
                    done += 1
                    _on_complete(done)
            else:
                with ProcessPoolExecutor(max_workers=args.workers) as pool:
                    futures = {
                        pool.submit(
                            process_one, job, str(per_gff_root), force=False
                        ): job
                        for job in prepared
                    }
                    for fut in as_completed(futures):
                        job = futures[fut]
                        try:
                            n_pairs, ann_id, skipped = fut.result()
                            if skipped:
                                skipped_existing += 1
                            else:
                                harvested += 1
                            if n_pairs == 0:
                                print(
                                    f"harvest empty {ann_id}: 0 xref pairs",
                                    file=sys.stderr,
                                )
                        except Exception as exc:  # noqa: BLE001
                            harvest_errors += 1
                            print(
                                f"harvest fail {job.annotation_id}: {exc}",
                                file=sys.stderr,
                            )
                        done += 1
                        _on_complete(done)

            print(
                f"harvested: ok={harvested} skipped_existing={skipped_existing} "
                f"errors={harvest_errors}",
                file=sys.stderr,
                flush=True,
            )

    # Attach any report id that is still absent from DB but has a ready shard
    # (includes harvest successes and pre-existing orphan shards).
    present_ids = set(known_ids)
    attach_candidates = sorted({r.annotation_id for r in new_rows})
    attached = 0
    attach_skipped = 0
    attach_errors = 0
    attach_missing_shard = 0
    for ann_id in attach_candidates:
        if ann_id in present_ids:
            attach_skipped += 1
            continue
        shard = shard_sqlite_path(per_gff_root, ann_id)
        if not shard.is_file():
            attach_missing_shard += 1
            print(
                f"attach skip: no shard for {ann_id}",
                file=sys.stderr,
            )
            continue
        try:
            key = next_annotation_key(conn)
            seed = read_seed_row_from_shard(per_gff_root, ann_id, key)
            ok = attach_shard(
                conn,
                per_gff_root,
                seed,
                taxonomy_parents=parents,
            )
            if ok:
                attached += 1
                present_ids.add(ann_id)
                print(
                    f"attach ok: {ann_id} → annotation_key={key}",
                    file=sys.stderr,
                    flush=True,
                )
            else:
                attach_skipped += 1
                present_ids.add(ann_id)
        except Exception as exc:  # noqa: BLE001
            attach_errors += 1
            print(f"attach fail {ann_id}: {exc}", file=sys.stderr)

    print(
        f"attach: ok={attached} skipped={attach_skipped} "
        f"missing_shard={attach_missing_shard} errors={attach_errors}",
        file=sys.stderr,
        flush=True,
    )

    refreshed, n_lineage = refresh_lineage_if_taxonomy_changed(
        conn, taxonomy_tsv, parents
    )
    if refreshed:
        print(
            f"lineage: refreshed ({n_lineage} rows) — taxonomy TSV changed",
            file=sys.stderr,
            flush=True,
        )
    else:
        print(
            f"lineage: unchanged ({n_lineage} rows) — taxonomy hash match",
            file=sys.stderr,
            flush=True,
        )

    n_stats = rebuild_namespace_stats(conn)
    conn.execute("ANALYZE")
    conn.commit()
    print(f"namespace_stats: {n_stats} namespace(s)", file=sys.stderr, flush=True)

    _touch_db_for_readers(db_path, conn)
    conn.close()

    if harvest_errors or attach_errors:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
