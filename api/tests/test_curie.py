"""Unit tests for CURIE parsing (sidecar §4)."""
from __future__ import annotations

import pytest

from helpers.curie import (
    CurieError,
    format_curie,
    normalize_go_accession,
    parse_curie,
)


def test_symbol_casefold() -> None:
    assert parse_curie("symbol:TP53") == ("symbol", "tp53")
    assert parse_curie("alias:Tp53") == ("alias", "tp53")


def test_go_padding() -> None:
    assert normalize_go_accession("8150") == "GO:0008150"
    assert parse_curie("GO:8150") == ("go", "GO:0008150")
    assert parse_curie("go:GO:0008150") == ("go", "GO:0008150")


def test_interpro_keeps_case() -> None:
    assert parse_curie("interpro:IPR002117") == ("interpro", "IPR002117")


def test_slash_namespace_short_prefix() -> None:
    assert parse_curie("swiss-prot:P04637") == ("uniprot/swiss-prot", "P04637")
    assert format_curie("uniprot/swiss-prot", "P04637") == "swiss-prot:P04637"


def test_unknown_prefix() -> None:
    with pytest.raises(CurieError, match="unknown CURIE prefix"):
        parse_curie("notans:foo")


def test_missing_colon() -> None:
    with pytest.raises(CurieError, match="colon"):
        parse_curie("tp53")
