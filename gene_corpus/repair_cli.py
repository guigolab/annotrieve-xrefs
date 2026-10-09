"""
Repair CLI: unwrap quote-wrapped symbol/alias accessions in gene_corpus.sqlite.

    # Dry-run (default): scan + write state files, no DB writes
    python -m gene_corpus.repair \\
        --work-dir /data/annotrieve/gene_corpus

    # Apply moves + recompute xref_meta / namespace_stats
    python -m gene_corpus.repair \\
        --work-dir /data/annotrieve/gene_corpus \\
        --apply

    # Recovery: recompute meta from affected_bare.txt after an interrupted apply
    python -m gene_corpus.repair \\
        --work-dir /data/annotrieve/gene_corpus \\
        --recompute-only

Copy/snapshot the DB before --apply on NFS. Prefer no concurrent writers.
API readers use immutable=1 and refresh on mtime. If apply dies after gene_hit
moves, re-run with --recompute-only.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.sql.repair_quotes import (
    DEFAULT_NAMESPACES,
    recompute_xref_meta_from_state,
    repair_wrapped_accessions,
)
from gene_corpus.sql.schema import connect_for_build


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Unwrap quote-wrapped symbol/alias accessions in "
            "gene_corpus.sqlite (gene_hit + xref_meta + namespace_stats). "
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
        "--apply",
        action="store_true",
        help="perform gene_hit moves and xref_meta recompute (default: dry-run)",
    )
    parser.add_argument(
        "--recompute-only",
        action="store_true",
        help=(
            "skip moves; recompute xref_meta from "
            "state-dir/affected_bare.txt (recovery after interrupted apply)"
        ),
    )
    parser.add_argument(
        "--namespaces",
        nargs="+",
        default=list(DEFAULT_NAMESPACES),
        metavar="NS",
        help=f"namespaces to repair (default: {' '.join(DEFAULT_NAMESPACES)})",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=None,
        help="directory for quoted.txt / affected_bare.txt "
        "(default: WORK_DIR/repair_quotes)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    db_path = work_dir / "gene_corpus.sqlite"
    state_dir: Path = (
        args.state_dir if args.state_dir is not None else work_dir / "repair_quotes"
    )
    apply: bool = bool(args.apply)
    recompute_only: bool = bool(args.recompute_only)
    namespaces: list[str] = list(args.namespaces)

    if apply and recompute_only:
        print(
            "--apply and --recompute-only cannot be combined",
            file=sys.stderr,
        )
        return 1
    if not namespaces:
        print("at least one --namespaces value is required", file=sys.stderr)
        return 1
    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1

    conn = connect_for_build(db_path)
    try:
        if recompute_only:
            result = recompute_xref_meta_from_state(conn, state_dir)
        else:
            result = repair_wrapped_accessions(
                conn,
                namespaces=namespaces,
                state_dir=state_dir,
                apply=apply,
            )
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        conn.close()
        return 1
    conn.close()

    mode = (
        "recompute-only"
        if recompute_only
        else ("apply" if result.applied else "dry-run")
    )
    print(
        f"repair {mode}: quoted_keys={result.quoted_keys} "
        f"affected_bare={result.affected_bare} "
        f"moved_rows={result.moved_rows} "
        f"dup_loci_skipped={result.dup_loci_skipped} "
        f"keys_repaired={result.keys_repaired} "
        f"wrapped_remaining={result.wrapped_remaining} "
        f"→ {db_path}",
        file=sys.stderr,
    )
    if result.applied and result.wrapped_remaining:
        print(
            f"repair warn: {result.wrapped_remaining} wrapped accession(s) "
            "still present",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
