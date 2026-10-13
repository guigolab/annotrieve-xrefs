"""
Index CLI: copy per-GFF tier_a_counts into xref_file.

    python -m gene_corpus.index --work-dir /data/gene_corpus
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.sql.schema import connect_for_build
from gene_corpus.sql.xref_index import index_missing_annotations


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build xref_file from per_gff/*/genes.sqlite tier_a_counts "
            "(SEARCHABLE_NAMESPACES only; presence for GET /annotations). "
            "Skips keys already in xref_file_status."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="corpus root (gene_corpus.sqlite + per_gff/)",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=50,
        metavar="N",
        help="WAL checkpoint every N newly indexed annotations (default: 50)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    per_gff_root = work_dir / "per_gff"
    db_path = work_dir / "gene_corpus.sqlite"

    if not db_path.is_file():
        print(
            f"gene_corpus.sqlite missing: {db_path} "
            "(run python -m gene_corpus.merge first)",
            file=sys.stderr,
        )
        return 1
    if not per_gff_root.is_dir():
        print(f"per_gff directory missing: {per_gff_root}", file=sys.stderr)
        return 1
    if args.checkpoint_every < 1:
        print("--checkpoint-every must be >= 1", file=sys.stderr)
        return 2

    conn = connect_for_build(db_path)
    try:
        placed, empty, errors = index_missing_annotations(
            conn,
            per_gff_root,
            checkpoint_every=args.checkpoint_every,
        )
        # TODO (reader mtime): when placed + empty > 0, mirror
        # sync_cli._touch_db_for_readers (WAL checkpoint + os.utime) so API
        # immutable=1 readers reopen after the index build.
    finally:
        conn.close()

    print(
        f"index done: placed={placed} empty={empty} errors={errors} "
        f"→ {db_path}",
        file=sys.stderr,
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
