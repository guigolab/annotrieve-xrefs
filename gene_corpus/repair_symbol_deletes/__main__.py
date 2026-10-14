"""Allow ``python -m gene_corpus.repair_symbol_deletes``."""
from __future__ import annotations

from gene_corpus.repair_symbol_deletes_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
