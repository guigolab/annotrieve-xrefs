"""
One-time CLI: add annotation.release_date + assembly_name and seed from
the annotations report.

    # Dry-run (default)
    python -m gene_corpus.seed_release_date --work-dir /data/gene_corpus

    # Apply ALTER + UPDATE
    python -m gene_corpus.seed_release_date --work-dir /data/gene_corpus --apply

    # Local TSV (must include release_date + assembly_name; skip download)
    python -m gene_corpus.seed_release_date --work-dir DIR --report PATH --apply
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.catalog.report import (
    DEFAULT_API_BASE,
    downloaded_annotation_report,
)
from gene_corpus.sql.schema import connect_for_build
from gene_corpus.sql.seed_release_date import (
    REPORT_SEED_FIELDS,
    parse_release_date_tsv,
    seed_release_dates,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "One-time job: ALTER annotation for release_date and "
            "assembly_name, then seed from POST /annotations/report "
            "(selected_fields=release_date; assembly_name is a default "
            "column). Default is dry-run; pass --apply to write."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="corpus root containing gene_corpus.sqlite",
    )
    parser.add_argument(
        "--api-base",
        default=DEFAULT_API_BASE,
        help=f"Annotrieve API base (default: {DEFAULT_API_BASE})",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help=(
            "local annotation report TSV with release_date and "
            "assembly_name (skip download)"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=600,
        help="timeout seconds for annotation report download",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help=(
            "ALTER table and UPDATE release_date + assembly_name "
            "(default: dry-run)"
        ),
    )
    return parser


def _run_seed(
    db_path: Path,
    report_path: Path,
    *,
    apply: bool,
):
    seeds, total, empty_date, empty_name = parse_release_date_tsv(report_path)
    print(
        f"report rows: {total}; with release_date: "
        f"{sum(1 for s in seeds.values() if s.release_date)}; "
        f"with assembly_name: "
        f"{sum(1 for s in seeds.values() if s.assembly_name)}; "
        f"empty release_date: {empty_date}; "
        f"empty assembly_name: {empty_name}",
        file=sys.stderr,
    )
    conn = connect_for_build(db_path)
    try:
        return seed_release_dates(
            conn,
            seeds,
            apply=apply,
            report_rows=total,
            empty_in_report=empty_date,
            empty_assembly_name=empty_name,
        )
    finally:
        conn.close()


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    db_path = work_dir / "gene_corpus.sqlite"
    apply: bool = bool(args.apply)

    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1

    if args.report is not None:
        if not args.report.is_file():
            print(f"report not found: {args.report}", file=sys.stderr)
            return 1
        print(f"report: {args.report}", file=sys.stderr)
        stats = _run_seed(db_path, args.report, apply=apply)
    else:
        print(
            f"downloading report from {args.api_base} "
            f"(selected_fields={REPORT_SEED_FIELDS})",
            file=sys.stderr,
        )
        with downloaded_annotation_report(
            args.api_base,
            selected_fields=REPORT_SEED_FIELDS,
            timeout=args.timeout,
        ) as report_path:
            stats = _run_seed(db_path, report_path, apply=apply)

    mode = "applied" if apply else "dry-run"
    added = ",".join(stats.columns_added) if stats.columns_added else "-"
    print(
        f"{mode}: columns_added={added}; "
        f"matched={stats.matched}; updated={stats.updated}; "
        f"report_ids_not_in_db={stats.missing_in_db}; "
        f"empty_release_date={stats.empty_release_date}; "
        f"empty_assembly_name={stats.empty_assembly_name}",
        file=sys.stderr,
    )
    if not apply:
        print("re-run with --apply to write", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
