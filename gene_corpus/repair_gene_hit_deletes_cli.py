"""
CLI: delete stale symbol/alias rows from gene_hit via xref_file anti-join.

    # Dry-run (count only)
    python -m gene_corpus.repair_gene_hit_deletes --corpus-db /data/gene_corpus.sqlite

    # Spill + delete
    python -m gene_corpus.repair_gene_hit_deletes --corpus-db /data/gene_corpus.sqlite --apply

    # Finalize only (scan done markers present)
    python -m gene_corpus.repair_gene_hit_deletes --corpus-db /data/gene_corpus.sqlite \\
        --apply --finalize-only
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.sql.repair_gene_hit_deletes import (
    NAMESPACES,
    XrefFileGuardError,
    cleanup_spill,
    default_spill_dir,
    finalize,
    is_scan_done,
    run_scan,
    scan_complete,
)
from gene_corpus.sql.schema import connect_for_build


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Delete gene_hit symbol/alias rows whose "
            "(namespace, accession, annotation_key) is missing from "
            "xref_file (annotation must be in xref_file_status). "
            f"Scan uses {len(NAMESPACES)} workers (one per namespace). "
            "Default is dry-run."
        )
    )
    parser.add_argument(
        "--corpus-db",
        type=Path,
        required=True,
        help="path to gene_corpus.sqlite",
    )
    parser.add_argument(
        "--spill-dir",
        type=Path,
        default=None,
        help=(
            "local TSV spill directory (default: "
            "annotrieve-xrefs/repair_gene_hit_deletes/)"
        ),
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write spill TSVs and delete from gene_hit (default: dry-run)",
    )
    parser.add_argument(
        "--finalize-only",
        action="store_true",
        help="skip scan; require part-{alias,symbol}.done",
    )
    parser.add_argument(
        "--keep-spill",
        action="store_true",
        help="do not delete spill dir after successful finalize",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="override empty xref_file guard",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    db_path: Path = args.corpus_db
    spill_dir: Path = (
        args.spill_dir if args.spill_dir is not None else default_spill_dir()
    )
    apply: bool = bool(args.apply)
    finalize_only: bool = bool(args.finalize_only)
    force: bool = bool(args.force)

    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1
    if finalize_only and not apply:
        print(
            "repair_gene_hit_deletes: --finalize-only requires --apply",
            file=sys.stderr,
        )
        return 1

    conn = connect_for_build(db_path)
    try:
        print(
            f"repair_gene_hit_deletes: spill_dir={spill_dir} "
            f"namespaces={','.join(NAMESPACES)} db={db_path}",
            file=sys.stderr,
            flush=True,
        )

        if not apply:
            try:
                summary = run_scan(
                    conn,
                    db_path=db_path,
                    spill_dir=spill_dir,
                    apply=False,
                    force=force,
                )
            except XrefFileGuardError as exc:
                print(f"repair_gene_hit_deletes: {exc}", file=sys.stderr)
                return 1
            for ns in summary.namespaces:
                samples = ", ".join(ns.sample_accessions[:10]) or "(none)"
                print(
                    f"repair_gene_hit_deletes dry-run {ns.namespace}: "
                    f"gene_hit={ns.gene_hit_rows:,} "
                    f"xref_file={ns.xref_file_rows:,} "
                    f"stale={ns.stale_rows:,} "
                    f"groups={ns.distinct_groups:,} "
                    f"samples=[{samples}]",
                    file=sys.stderr,
                )
            print(
                f"repair_gene_hit_deletes dry-run total stale="
                f"{summary.stale_rows:,} → {db_path}",
                file=sys.stderr,
            )
            print(
                "pass --apply to spill deletes and remove gene_hit rows",
                file=sys.stderr,
            )
            return 0

        if not finalize_only:
            try:
                summary = run_scan(
                    conn,
                    db_path=db_path,
                    spill_dir=spill_dir,
                    apply=True,
                    force=force,
                )
            except XrefFileGuardError as exc:
                print(f"repair_gene_hit_deletes: {exc}", file=sys.stderr)
                return 1
            for ns in summary.namespaces:
                print(
                    f"repair_gene_hit_deletes scan {ns.namespace}: "
                    f"stale={ns.stale_rows:,} groups={ns.distinct_groups:,}",
                    file=sys.stderr,
                    flush=True,
                )
            print(
                f"repair_gene_hit_deletes scan: "
                f"tsv_lines={summary.tsv_lines_appended:,}",
                file=sys.stderr,
                flush=True,
            )
        else:
            if not scan_complete(spill_dir):
                done = {
                    ns: is_scan_done(spill_dir, ns) for ns in NAMESPACES
                }
                print(
                    "repair_gene_hit_deletes: --finalize-only but scan "
                    f"incomplete (done={done})",
                    file=sys.stderr,
                )
                return 1

        if not scan_complete(spill_dir):
            done = {ns: is_scan_done(spill_dir, ns) for ns in NAMESPACES}
            print(
                f"repair_gene_hit_deletes: scan incomplete (done={done}); "
                "re-run --apply to resume",
                file=sys.stderr,
            )
            return 0

        fin = finalize(conn, spill_dir=spill_dir)
        print(
            f"repair_gene_hit_deletes finalize: "
            f"tsv_lines={fin.tsv_lines:,} deleted={fin.deleted:,} "
            f"skipped_resume={fin.skipped_resume:,}",
            file=sys.stderr,
            flush=True,
        )
        if not args.keep_spill:
            cleanup_spill(spill_dir)
            print(
                f"repair_gene_hit_deletes: cleaned {spill_dir}",
                file=sys.stderr,
            )
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
