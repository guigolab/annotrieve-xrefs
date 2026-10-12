"""Unit tests for source_database → integer codes."""
from __future__ import annotations

from helpers.source_codes import source_code


def test_source_code_known() -> None:
    assert source_code("Ensembl") == 1
    assert source_code("ensembl") == 1
    assert source_code("GenBank") == 2
    assert source_code("RefSeq") == 3
    assert source_code("refseq") == 3


def test_source_code_unknown() -> None:
    assert source_code(None) == 0
    assert source_code("") == 0
    assert source_code("Other") == 0
