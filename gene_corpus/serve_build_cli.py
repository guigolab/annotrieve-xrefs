"""
Build the served presence index from the global gene_corpus DB.

    python -m gene_corpus.serve_build --work-dir /data/gene_corpus
    python -m gene_corpus.serve_build --work-dir /data/gene_corpus --name-prefix tp5 --out-dir /tmp/trial

Writes ``<out-dir>/v<data_version>/presence.sqlite`` and switches the
``<out-dir>/current`` symlink atomically. Partial builds (``--name-prefix`` /
``--limit-names``) are written to ``<out-dir>/partial-<data_version>/`` and never
switch ``current``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from gene_corpus.sql.serve_presence import BuildError, build_presence, publish_version


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build presence.sqlite (compact gene-name index) from gene_corpus.sqlite."
    )
    parser.add_argument("--work-dir", type=Path, required=True,
                        help="corpus root containing gene_corpus.sqlite")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="default: <work-dir>/serve")
    parser.add_argument("--name-prefix", default=None,
                        help="trial build: only names starting with this prefix")
    parser.add_argument("--limit-names", type=int, default=None,
                        help="trial build: stop after N names")
    parser.add_argument("--verify-sample", type=int, default=200,
                        help="names re-checked against xref_file after the build (0 = skip)")
    parser.add_argument("--tmp-dir", type=Path, default=None,
                        help="scratch dir for SQLite temp files")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    src = args.work_dir / "gene_corpus.sqlite"
    if not src.is_file():
        print(f"gene_corpus.sqlite missing: {src}", file=sys.stderr)
        return 2
    out_dir = args.out_dir or (args.work_dir / "serve")
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.tmp_dir:
        os.environ["SQLITE_TMPDIR"] = str(args.tmp_dir)

    building = out_dir / ".presence.building.sqlite"
    try:
        stats = build_presence(
            src,
            building,
            name_prefix=args.name_prefix,
            limit_names=args.limit_names,
            verify_sample=args.verify_sample,
        )
    except BuildError as exc:
        print(f"BUILD FAILED: {exc}", file=sys.stderr)
        if building.exists():
            building.unlink()
        return 1

    if stats["partial"]:
        part = out_dir / f"partial-{stats['data_version']}"
        part.mkdir(parents=True, exist_ok=True)
        final = part / "presence.sqlite"
        os.replace(building, final)
        (part / "build_report.json").write_text(json.dumps(stats, indent=2))
        print(f"partial build written to {final} (current not changed)", file=sys.stderr)
    else:
        final = publish_version(out_dir, building, stats)
        print(f"published {final}; current -> {(out_dir / 'current').resolve().name}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
