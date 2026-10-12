"""
Root CLI dispatcher — use the split harvest / merge / sync / repair / labels
entry points.

    python -m gene_corpus.harvest --help
    python -m gene_corpus.merge --help
    python -m gene_corpus.sync --help
    python -m gene_corpus.repair --help
    python -m gene_corpus.labels --help
"""
from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    _ = argv
    print(
        "gene_corpus: use separate harvest, merge, sync, repair, and labels "
        "commands:\n"
        "  python -m gene_corpus.harvest --work-dir DIR "
        "--files-root FILES [--annotation-id ID ...]\n"
        "  python -m gene_corpus.merge --work-dir DIR "
        "--taxonomy-tsv PATH\n"
        "  python -m gene_corpus.sync --work-dir DIR "
        "--files-root FILES [--taxonomy-tsv PATH]\n"
        "  python -m gene_corpus.repair --work-dir DIR [--apply]\n"
        "  python -m gene_corpus.labels --work-dir DIR [--namespace go]\n",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
