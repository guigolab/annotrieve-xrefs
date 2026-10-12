"""Registered label providers (namespace → provider)."""
from __future__ import annotations

from gene_corpus.labels.providers.base import LabelProvider, LabelTerm
from gene_corpus.labels.providers.go import GoLabelProvider
from gene_corpus.labels.providers.interpro import InterproLabelProvider
from gene_corpus.labels.providers.mirbase import MirbaseLabelProvider
from gene_corpus.labels.providers.pfam import PfamLabelProvider
from gene_corpus.labels.providers.rfam import RfamLabelProvider
from gene_corpus.labels.providers.swissprot import SwissprotLabelProvider

PROVIDERS: dict[str, LabelProvider] = {
    "go": GoLabelProvider(),
    "interpro": InterproLabelProvider(),
    "mirbase": MirbaseLabelProvider(),
    "pfam": PfamLabelProvider(),
    "rfam": RfamLabelProvider(),
    "uniprot/swiss-prot": SwissprotLabelProvider(),
}

__all__ = [
    "PROVIDERS",
    "LabelProvider",
    "LabelTerm",
    "GoLabelProvider",
    "InterproLabelProvider",
    "MirbaseLabelProvider",
    "PfamLabelProvider",
    "RfamLabelProvider",
    "SwissprotLabelProvider",
]
