"""
CLI: delete-delta cleanup of via_level>0 symbol/alias.

    # Dry-run
    python -m gene_corpus.repair_symbol_deletes --work-dir /data/gene_corpus

    # Spill deletes to local part TSVs, then patch xref_file / xref_meta
    python -m gene_corpus.repair_symbol_deletes --work-dir /data/gene_corpus --apply

    # Finalize only (scan cursors complete)
    python -m gene_corpus.repair_symbol_deletes --work-dir /data/gene_corpus \\
        --apply --finalize-only
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.sql.repair_symbol_deletes import (
    N_WORKERS,
    cleanup_spill,
    default_spill_dir,
    finalize,
    list_annotations,
    read_cursor,
    run_scan,
    scan_complete,
)
from gene_corpus.sql.schema import connect_for_build


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Delete via_level>0 symbol/alias from per-GFF shards, spill "
            "delete rows to local part-0/1.tsv under the annotrieve-xrefs "
            "tree, then patch xref_file / xref_meta. Always uses "
            f"{N_WORKERS} workers. Default is dry-run."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="corpus root (gene_corpus.sqlite + per_gff/)",
    )
    parser.add_argument(
        "--spill-dir",
        type=Path,
        default=None,
        help=(
            "local TSV spill directory (default: "
            "annotrieve-xrefs/repair_symbol_deletes/)"
        ),
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write shards, part TSVs, and finalize (default: dry-run)",
    )
    parser.add_argument(
        "--finalize-only",
        action="store_true",
        help="skip scan; require partition cursors complete",
    )
    parser.add_argument(
        "--keep-spill",
        action="store_true",
        help="do not delete spill dir after successful finalize",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    spill_dir: Path = (
        args.spill_dir if args.spill_dir is not None else default_spill_dir()
    )
    db_path = work_dir / "gene_corpus.sqlite"
    per_gff_root = work_dir / "per_gff"
    apply: bool = bool(args.apply)
    finalize_only: bool = bool(args.finalize_only)

    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1
    if not per_gff_root.is_dir():
        print(f"per_gff directory missing: {per_gff_root}", file=sys.stderr)
        return 1

    conn = connect_for_build(db_path)
    try:
        annotations = list_annotations(conn)
        print(
            f"repair_symbol_deletes: spill_dir={spill_dir} "
            f"workers={N_WORKERS} annotations={len(annotations)}",
            file=sys.stderr,
            flush=True,
        )

        if not apply:
            summary = run_scan(
                conn,
                work_dir=work_dir,
                spill_dir=spill_dir,
                apply=False,
            )
            print(
                f"repair_symbol_deletes dry-run: "
                f"annotations={summary.total_annotations} "
                f"would_delete_xref_rows≈{summary.xref_rows_deleted} "
                f"→ {db_path}",
                file=sys.stderr,
            )
            print(
                "pass --apply to spill deletes, clean shards, and finalize",
                file=sys.stderr,
            )
            return 0

        if not finalize_only:
            summary = run_scan(
                conn,
                work_dir=work_dir,
                spill_dir=spill_dir,
                apply=True,
            )
            print(
                f"repair_symbol_deletes scan: "
                f"processed={summary.processed} "
                f"skipped_resume={summary.skipped_resume} "
                f"missing_shard={summary.missing_shard} "
                f"deleted_xref={summary.xref_rows_deleted} "
                f"tsv_lines={summary.tsv_lines_appended} "
                f"total={summary.total_annotations}",
                file=sys.stderr,
                flush=True,
            )
        else:
            if not scan_complete(conn, spill_dir):
                cursors = [
                    read_cursor(spill_dir, w) for w in range(N_WORKERS)
                ]
                print(
                    "repair_symbol_deletes: --finalize-only but scan "
                    f"incomplete (cursors={cursors} "
                    f"total={len(annotations)})",
                    file=sys.stderr,
                )
                return 1

        if not scan_complete(conn, spill_dir):
            cursors = [read_cursor(spill_dir, w) for w in range(N_WORKERS)]
            print(
                f"repair_symbol_deletes: scan incomplete "
                f"(cursors={cursors} total={len(annotations)}); "
                "re-run --apply to resume",
                file=sys.stderr,
            )
            return 0

        fin = finalize(conn, spill_dir=spill_dir)
        print(
            f"repair_symbol_deletes finalize: "
            f"tsv_lines={fin.tsv_lines} "
            f"unique_deletes={fin.unique_deletes} "
            f"xref_file_deleted={fin.xref_file_deleted} "
            f"xref_meta_updated={fin.xref_meta_updated} "
            f"xref_meta_deleted={fin.xref_meta_deleted}",
            file=sys.stderr,
            flush=True,
        )
        if not args.keep_spill:
            cleanup_spill(spill_dir)
            print(
                f"repair_symbol_deletes: cleaned {spill_dir}",
                file=sys.stderr,
            )
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
