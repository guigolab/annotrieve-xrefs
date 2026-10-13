"""Allow ``python -m gene_corpus.prune_xref``."""
from __future__ import annotations

from gene_corpus.prune_xref_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
