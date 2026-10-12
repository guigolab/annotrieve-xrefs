"""
Labels CLI: fill xref_label / xref_label_alias and refresh xref_label_fts.

    python -m gene_corpus.labels \\
        --work-dir /data/gene_corpus \\
        --namespace go

Each load also rebuilds the trigram FTS index (``xref_label_fts``) for that
namespace so labels/aliases are searchable by substring.

Re-run after sync when new xref_meta accessions appear (sync wiring later).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers import PROVIDERS


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Load human-readable labels for xref_meta accessions "
            "(GO, miRBase, Rfam, InterPro, Pfam, UniProt/Swiss-Prot) "
            "and refresh xref_label_fts."
        )
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="corpus root containing gene_corpus.sqlite",
    )
    parser.add_argument(
        "--namespace",
        action="append",
        dest="namespaces",
        default=None,
        metavar="NS",
        help=(
            "namespace to load (repeatable; default: all registered "
            f"providers: {', '.join(sorted(PROVIDERS))})"
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="ontology download cache (default: WORK_DIR/label_cache)",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-download ontology files even if cached",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    work_dir: Path = args.work_dir
    namespaces: list[str] = list(args.namespaces or sorted(PROVIDERS))

    unknown = [ns for ns in namespaces if ns not in PROVIDERS]
    if unknown:
        print(
            f"unknown label namespace(s): {', '.join(unknown)}; "
            f"registered: {', '.join(sorted(PROVIDERS))}",
            file=sys.stderr,
        )
        return 1

    db_path = work_dir / "gene_corpus.sqlite"
    if not db_path.is_file():
        print(f"gene_corpus.sqlite missing: {db_path}", file=sys.stderr)
        return 1

    for namespace in namespaces:
        try:
            load_namespace_labels(
                work_dir,
                namespace,
                cache_dir=args.cache_dir,
                force_download=bool(args.force_download),
            )
        except Exception as exc:
            print(f"labels {namespace}: failed: {exc}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
