"""Allow ``python -m gene_corpus.labels``."""
from __future__ import annotations

from gene_corpus.labels.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
