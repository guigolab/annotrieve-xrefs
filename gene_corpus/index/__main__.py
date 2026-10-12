"""Allow ``python -m gene_corpus.index``."""
from __future__ import annotations

from gene_corpus.index_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
