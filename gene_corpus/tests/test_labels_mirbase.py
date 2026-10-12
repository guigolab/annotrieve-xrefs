"""miRBase label parse + RNAcentral fallback + load against xref_meta."""
from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.mirbase import (
    MirbaseLabelProvider,
    _iter_official,
    _iter_rnacentral_tsv,
)
from gene_corpus.sql.schema import connect_for_build, init_schema


def _write_official_fixture(cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "mirna.txt").write_text(
        "1\tMI0001445\thsa-mir-12136\t\tdesc\n"
        "2\tMI0006387\thsa-mir-6859-1\t\tdesc\n",
        encoding="utf-8",
    )
    (cache / "mirna_mature.txt").write_text(
        "10\thsa-miR-12136\tMIMAT0031928\t\n"
        "11\thsa-miR-6859-5p\tMIMAT0026718\t\n",
        encoding="utf-8",
    )
    (cache / "dead_mirna.txt").write_text(
        "MI0009656\told-mir-dead\t\n",
        encoding="utf-8",
    )
    (cache / "aliases.txt").write_text(
        "MI0001445\thsa-mir-12136;hsa-mir-12136-old\n"
        "MIMAT0031928\thsa-miR-12136;hsa-miR-12136-alt\n",
        encoding="utf-8",
    )
    (cache / "source.txt").write_text("official\n", encoding="utf-8")


def _write_rnacentral_fixture(path: Path) -> None:
    path.write_text(
        "URS1\tMIRBASE\tMI0001445\t9606\tpre_miRNA\thsa-mir-12136\n"
        "URS2\tMIRBASE\tMIMAT0031928\t9606\tmiRNA\thsa-miR-12136\n"
        "URS3\tMIRBASE\tMI0006387\t9606\tpre_miRNA\thsa-mir-6859-1\n",
        encoding="utf-8",
    )


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('mirbase', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_iter_official_maps_mi_mimat_dead_and_aliases(tmp_path: Path) -> None:
    cache = tmp_path / "mirbase"
    _write_official_fixture(cache)
    by_id = {t.accession: t for t in _iter_official(cache)}
    assert by_id["MI0001445"].label == "hsa-mir-12136"
    assert "hsa-mir-12136-old" in by_id["MI0001445"].aliases
    assert by_id["MIMAT0031928"].label == "hsa-miR-12136"
    assert "hsa-miR-12136-alt" in by_id["MIMAT0031928"].aliases
    assert by_id["MI0009656"].label == "old-mir-dead"


def test_iter_rnacentral_fallback(tmp_path: Path) -> None:
    tsv = tmp_path / "rnacentral_mirbase.tsv"
    _write_rnacentral_fixture(tsv)
    terms = list(_iter_rnacentral_tsv(tsv))
    by_id = {t.accession: t for t in terms}
    assert by_id["MI0001445"].label == "hsa-mir-12136"
    assert by_id["MIMAT0031928"].label == "hsa-miR-12136"
    assert by_id["MI0001445"].aliases == ()


def test_iter_terms_uses_source_marker(tmp_path: Path) -> None:
    cache = tmp_path / "mirbase"
    cache.mkdir()
    _write_rnacentral_fixture(cache / "rnacentral_mirbase.tsv")
    (cache / "source.txt").write_text("rnacentral\n", encoding="utf-8")
    terms = list(MirbaseLabelProvider().iter_terms(cache))
    assert {t.accession for t in terms} == {
        "MI0001445",
        "MIMAT0031928",
        "MI0006387",
    }


def test_download_falls_back_to_rnacentral(tmp_path: Path, monkeypatch) -> None:
    provider = MirbaseLabelProvider()
    cache_root = tmp_path / "label_cache"

    def boom(url: str, dest: Path) -> None:
        if "rnacentral" in url or "RNAcentral" in url or "mirbase.tsv" in url:
            _write_rnacentral_fixture(dest)
            return
        raise OSError("mirbase.org down")

    monkeypatch.setattr(
        "gene_corpus.labels.providers.mirbase._download_to", boom
    )
    path = provider.download(cache_root, force=True)
    assert (path / "source.txt").read_text(encoding="utf-8").strip() == (
        "rnacentral"
    )
    terms = list(provider.iter_terms(path))
    assert any(t.accession == "MI0001445" for t in terms)


def test_load_writes_labels_from_rnacentral_fixture() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["MI0001445", "MIMAT0031928", "MI9999999"])
        conn.close()

        src = work / "rnacentral_mirbase.tsv"
        _write_rnacentral_fixture(src)
        stats = load_namespace_labels(
            work, "mirbase", source_path=src, force_download=False
        )
        assert stats.labels_written == 2
        assert stats.meta_missing_label == 1

        conn = sqlite3.connect(str(db))
        labels = {
            (acc, lab)
            for acc, lab in conn.execute(
                "SELECT accession, label FROM xref_label "
                "WHERE namespace='mirbase'"
            )
        }
        assert ("MI0001445", "hsa-mir-12136") in labels
        assert ("MIMAT0031928", "hsa-miR-12136") in labels
        conn.close()
