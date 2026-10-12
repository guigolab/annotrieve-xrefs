"""Rfam label parse + load against xref_meta."""
from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.rfam import (
    RfamLabelProvider,
    _iter_dead_file,
    _iter_family_file,
)
from gene_corpus.sql.schema import connect_for_build, init_schema


def _family_row(
    acc: str,
    rfam_id: str,
    description: str,
    previous_id: str = "NULL",
) -> str:
    # Pad to col 10 (previous_id) matching family.txt dump layout.
    cells = [""] * 11
    cells[0] = acc
    cells[1] = rfam_id
    cells[2] = "1"
    cells[3] = description
    cells[10] = previous_id
    return "\t".join(cells) + "\n"


def _write_rfam_fixture(cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "family.txt").write_text(
        _family_row("RF00001", "5S_rRNA", "5S ribosomal RNA")
        + _family_row(
            "RF00005",
            "tRNA",
            "tRNA",
            previous_id="tRNA_old;trna",
        ),
        encoding="utf-8",
    )
    (cache / "dead_family.txt").write_text(
        "RF00098\tunknown\tsuspected pseudogene\tNULL\tDead title\tjgt\n",
        encoding="utf-8",
    )


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('rfam', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_iter_family_maps_id_description_previous(tmp_path: Path) -> None:
    path = tmp_path / "family.txt"
    path.write_text(
        _family_row("RF00001", "5S_rRNA", "5S ribosomal RNA")
        + _family_row("RF00005", "tRNA", "tRNA", previous_id="tRNA_old;trna"),
        encoding="utf-8",
    )
    by_id = {t.accession: t for t in _iter_family_file(path)}
    assert by_id["RF00001"].label == "5S_rRNA"
    assert by_id["RF00001"].aliases == ("5S ribosomal RNA",)
    assert by_id["RF00005"].label == "tRNA"
    # description equals label → skipped; previous_id tokens kept
    assert "tRNA_old" in by_id["RF00005"].aliases
    assert "trna" in by_id["RF00005"].aliases


def test_iter_dead_keeps_title_ignores_forward(tmp_path: Path) -> None:
    path = tmp_path / "dead_family.txt"
    path.write_text(
        "RF00098\tunknown\tcomment\tRF00001\tDead title\tjgt\n",
        encoding="utf-8",
    )
    terms = list(_iter_dead_file(path))
    assert len(terms) == 1
    assert terms[0].accession == "RF00098"
    assert terms[0].label == "unknown"
    assert terms[0].aliases == ("Dead title",)
    assert terms[0].alt_ids == ()


def test_iter_terms_directory(tmp_path: Path) -> None:
    cache = tmp_path / "rfam"
    _write_rfam_fixture(cache)
    by_id = {t.accession: t for t in RfamLabelProvider().iter_terms(cache)}
    assert "RF00001" in by_id
    assert "RF00098" in by_id


def test_load_writes_labels_for_meta_ids() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["RF00001", "RF00098", "RF99999"])
        conn.close()

        cache = work / "fixture_rfam"
        _write_rfam_fixture(cache)
        stats = load_namespace_labels(
            work, "rfam", source_path=cache, force_download=False
        )
        assert stats.labels_written == 2
        assert stats.meta_missing_label == 1

        conn = sqlite3.connect(str(db))
        labels = {
            (acc, lab)
            for acc, lab in conn.execute(
                "SELECT accession, label FROM xref_label WHERE namespace='rfam'"
            )
        }
        assert ("RF00001", "5S_rRNA") in labels
        assert ("RF00098", "unknown") in labels
        aliases = {
            alias
            for (alias,) in conn.execute(
                "SELECT alias FROM xref_label_alias "
                "WHERE namespace='rfam' AND accession='RF00001'"
            )
        }
        assert "5S ribosomal RNA" in aliases
        conn.close()
