"""
One-time CLI: add annotation.release_date and seed from the annotations report.

    # Dry-run (default)
    python -m gene_corpus.seed_release_date --work-dir /data/gene_corpus

    # Apply ALTER + UPDATE
    python -m gene_corpus.seed_release_date --work-dir /data/gene_corpus --apply

    # Local TSV (must include release_date; skip download)
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
    REPORT_RELEASE_DATE_FIELDS,
    parse_release_date_tsv,
    seed_release_dates,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "One-time job: ALTER annotation for release_date and seed values "
            "from POST /annotations/report?selected_fields=release_date. "
            "Default is dry-run; pass --apply to write."
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
        help="local annotation report TSV with release_date (skip download)",
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
        help="ALTER table and UPDATE release_date (default: dry-run)",
    )
    return parser


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
        id_to_date, total, empty = parse_release_date_tsv(args.report)
        print(
            f"report: {args.report} ({total} rows, "
            f"{len(id_to_date)} with release_date, {empty} empty)",
            file=sys.stderr,
        )
        conn = connect_for_build(db_path)
        try:
            stats = seed_release_dates(
                conn,
                id_to_date,
                apply=apply,
                report_rows=total,
                empty_in_report=empty,
            )
        finally:
            conn.close()
    else:
        print(
            f"downloading report from {args.api_base} "
            f"(selected_fields={REPORT_RELEASE_DATE_FIELDS})",
            file=sys.stderr,
        )
        with downloaded_annotation_report(
            args.api_base,
            selected_fields=REPORT_RELEASE_DATE_FIELDS,
            timeout=args.timeout,
        ) as report_path:
            id_to_date, total, empty = parse_release_date_tsv(report_path)
            print(
                f"report rows: {total}; with release_date: {len(id_to_date)}; "
                f"empty: {empty}",
                file=sys.stderr,
            )
            conn = connect_for_build(db_path)
            try:
                stats = seed_release_dates(
                    conn,
                    id_to_date,
                    apply=apply,
                    report_rows=total,
                    empty_in_report=empty,
                )
            finally:
                conn.close()

    mode = "applied" if apply else "dry-run"
    print(
        f"{mode}: column_added={stats.column_added}; "
        f"matched={stats.matched}; updated={stats.updated}; "
        f"report_ids_not_in_db={stats.missing_in_db}; "
        f"empty_in_report={stats.empty_in_report}",
        file=sys.stderr,
    )
    if not apply:
        print("re-run with --apply to write", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
