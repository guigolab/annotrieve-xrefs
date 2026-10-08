"""Unit tests for GFF strand formatting."""
from __future__ import annotations

from helpers.strand import format_strand


def test_format_strand() -> None:
    assert format_strand(1) == "+"
    assert format_strand(-1) == "-"
    assert format_strand(0) == "."
    assert format_strand(99) == "."
