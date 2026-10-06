"""Allow ``python -m gene_corpus.sync``."""
from __future__ import annotations

from gene_corpus.sync_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
