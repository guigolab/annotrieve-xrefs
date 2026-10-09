"""Allow ``python -m gene_corpus.repair``."""
from __future__ import annotations

from gene_corpus.repair_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
