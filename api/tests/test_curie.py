"""Unit tests for CURIE parsing (sidecar §4)."""
from __future__ import annotations

import pytest

from helpers.curie import (
    CurieError,
    normalize_go_accession,
    parse_curie_local,
    resolve_namespace_local,
)


def test_symbol_casefold() -> None:
    assert parse_curie_local("symbol:TP53") == ("symbol", "tp53")
    assert parse_curie_local("alias:Tp53") == ("alias", "tp53")


def test_go_padding() -> None:
    assert normalize_go_accession("8150") == "GO:0008150"
    assert parse_curie_local("GO:8150") == ("go", "GO:0008150")
    assert parse_curie_local("go:GO:0008150") == ("go", "GO:0008150")


def test_interpro_keeps_case() -> None:
    assert parse_curie_local("interpro:IPR002117") == ("interpro", "IPR002117")


def test_slash_namespace_short_prefix() -> None:
    assert parse_curie_local("swiss-prot:P04637") == ("uniprot/swiss-prot", "P04637")


def test_missing_colon() -> None:
    with pytest.raises(CurieError, match="colon"):
        parse_curie_local("tp53")


def test_parse_curie_local_tier_b() -> None:
    assert parse_curie_local("refseq_transcript:XM_015139305.2") == (
        "refseq_transcript",
        "XM_015139305.2",
    )


def test_parse_curie_local_swiss_prot() -> None:
    assert parse_curie_local("swiss-prot:P04637") == (
        "uniprot/swiss-prot",
        "P04637",
    )


def test_resolve_namespace_local_slash_form() -> None:
    assert resolve_namespace_local("uniprot/swiss-prot") == "uniprot/swiss-prot"
    assert resolve_namespace_local("swiss-prot") == "uniprot/swiss-prot"
    assert resolve_namespace_local("") is None
