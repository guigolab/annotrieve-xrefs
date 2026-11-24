"""
gene_corpus: per-annotation gene shards and the served search index built from them.

    python -m gene_corpus.harvest      # annotation report -> per-annotation genes.sqlite shards
    python -m gene_corpus.serve_build  # shards -> sidecars -> presence / loci / loci_ref files
"""

from gene_corpus.namespaces import canonicalize_db, include_in_local, include_in_v1
from gene_corpus.profiles import PROFILES

__all__ = [
    "PROFILES",
    "canonicalize_db",
    "include_in_local",
    "include_in_v1",
]
