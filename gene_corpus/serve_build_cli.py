"""
Build the served index files from the per-annotation shards.

    python -m gene_corpus.serve_build catalogue  --work-dir DIR   # annotation metadata -> manifest
    python -m gene_corpus.serve_build sidecars   --work-dir DIR   # one sidecar per shard
    python -m gene_corpus.serve_build build      --work-dir DIR   # presence.sqlite, published as current
    python -m gene_corpus.serve_build loci       --work-dir DIR   # loci.sqlite next to it
    python -m gene_corpus.serve_build loci-ref   --work-dir DIR   # loci_ref.sqlite (locus_ref lookup)
    python -m gene_corpus.serve_build missing|status --work-dir DIR

See ``serve_backfill_cli`` for the options of each sub-command.
"""
from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    from gene_corpus.serve_backfill_cli import main as backfill_main

    return backfill_main(list(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    raise SystemExit(main())
