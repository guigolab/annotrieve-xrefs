"""
Prune CLI: keep only SEARCHABLE_NAMESPACES in xref_file, then VACUUM.

    # Dry-run (default)
    python -m gene_corpus.prune_xref --work-dir /data/gene_corpus

    # Delete demoted rows, fix placed, reclaim disk
    python -m gene_corpus.prune_xref --work-dir /data/gene_corpus --apply
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.namespaces.tiers import SEARCHABLE_NAMESPACES
from gene_corpus.sql.prune_xref_file import apply_prune, summarize_prune
from gene_corpus.sql.schema import connect_for_build


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Delete xref_file rows outside SEARCHABLE_NAMESPACES, recompute "
            "xref_file_status.placed, and VACUUM. Default is dry-run; pass "
            "--apply to write. Prefer no concurrent writers."
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
        help="delete demoted rows, fix placed, and VACUUM (default: dry-run)",
    )
    return parser


def _fmt_bytes(n: int) -> str:
    if n < 1024:
        return f"{n} B"
    if n < 1024**2:
        return f"{n / 1024:.1f} KiB"
    if n < 1024**3:
        return f"{n / 1024**2:.1f} MiB"
    return f"{n / 1024**3:.2f} GiB"


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    db_path = work_dir / "gene_corpus.sqlite"
    apply: bool = bool(args.apply)

    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1

    conn = connect_for_build(db_path)
    try:
        if apply:
            result = apply_prune(conn, db_path=db_path)
            summary = result.summary
            print(
                f"prune_xref apply: deleted={result.rows_deleted} "
                f"kept={summary.rows_kept} "
                f"placed_flipped={result.placed_flipped} "
                f"vacuum={result.vacuumed} "
                f"→ {db_path}",
                file=sys.stderr,
            )
            if result.size_before is not None and result.size_after is not None:
                print(
                    f"prune_xref size: {_fmt_bytes(result.size_before)} → "
                    f"{_fmt_bytes(result.size_after)}",
                    file=sys.stderr,
                )
        else:
            summary = summarize_prune(conn)
            print(
                f"prune_xref dry-run: rows_total={summary.rows_total} "
                f"to_delete={summary.rows_to_delete} "
                f"kept={summary.rows_kept} "
                f"placed_to_zero={summary.placed_to_zero} "
                f"(allowlist={len(SEARCHABLE_NAMESPACES)} ns) "
                f"→ {db_path}",
                file=sys.stderr,
            )
            if summary.namespaces_to_delete:
                for ns, n in summary.namespaces_to_delete[:20]:
                    print(f"  {ns}: {n}", file=sys.stderr)
                if len(summary.namespaces_to_delete) > 20:
                    print(
                        f"  … {len(summary.namespaces_to_delete) - 20} more",
                        file=sys.stderr,
                    )
            print(
                "pass --apply to delete, fix placed, and VACUUM",
                file=sys.stderr,
            )
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
