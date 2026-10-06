"""Unit tests for annotation_lineage helpers."""
from __future__ import annotations

from pathlib import Path

import pytest

from gene_corpus.sql.lineage import lineage_taxids, load_parent_map


def test_load_parent_map(tmp_path: Path) -> None:
    tsv = tmp_path / "t.tsv"
    tsv.write_text(
        "taxid\tparent_taxid\tscientific_name\n"
        "9606\t9605\tHomo sapiens\n"
        "9605\t\tHomo\n",
        encoding="utf-8",
    )
    parents = load_parent_map(tsv)
    assert parents[9606] == 9605
    assert parents[9605] is None


def test_load_parent_map_requires_columns(tmp_path: Path) -> None:
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="taxid and parent_taxid"):
        load_parent_map(tsv)


def test_lineage_stops_on_cycle() -> None:
    parents = {1: 2, 2: 1}
    assert lineage_taxids(1, parents) == [1, 2]
