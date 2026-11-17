"""Allow ``python -m gene_corpus.serve_build``."""
from __future__ import annotations

from gene_corpus.serve_build_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
