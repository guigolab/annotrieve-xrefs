"""Allow ``python -m gene_corpus.seed_release_date``."""
from __future__ import annotations

from gene_corpus.seed_release_date_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
