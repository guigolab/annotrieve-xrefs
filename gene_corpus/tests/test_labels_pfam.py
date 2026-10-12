"""Pfam label parse + load against xref_meta."""
from __future__ import annotations

import gzip
import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.pfam import (
    PfamLabelProvider,
    _iter_clans_tsv,
    normalize_pfam_accession,
)
from gene_corpus.sql.schema import connect_for_build, init_schema

_CLANS_FIXTURE = """\
PF00001	CL0192	GPCR_A	7tm_1	7 transmembrane receptor (rhodopsin family)
PF00008	CL0001	EGF	EGF	EGF-like domain
PF00069	CL0016	PKinase	Pkinase	Protein kinase domain
"""


def _write_clans(path: Path, *, gzipped: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = _CLANS_FIXTURE.encode("utf-8")
    if gzipped:
        path.write_bytes(gzip.compress(data))
    else:
        path.write_bytes(data)


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('pfam', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_normalize_pfam_strips_version() -> None:
    assert normalize_pfam_accession("PF00001") == "PF00001"
    assert normalize_pfam_accession("PF00001.15") == "PF00001"
    assert normalize_pfam_accession("pf00069.1") == "PF00069"
    assert normalize_pfam_accession("IPR000008,InterPro") is None


def test_iter_clans_maps_description_and_short_id(tmp_path: Path) -> None:
    path = tmp_path / "Pfam-A.clans.tsv"
    _write_clans(path)
    by_id = {t.accession: t for t in _iter_clans_tsv(path)}
    assert by_id["PF00001"].label == (
        "7 transmembrane receptor (rhodopsin family)"
    )
    assert by_id["PF00001"].aliases == ("7tm_1",)
    assert "GPCR_A" not in by_id["PF00001"].aliases
    assert "CL0192" not in by_id["PF00001"].aliases
    # short id equals a token in label path: EGF == not full label
    assert by_id["PF00008"].label == "EGF-like domain"
    assert by_id["PF00008"].aliases == ("EGF",)


def test_iter_terms_gzipped(tmp_path: Path) -> None:
    path = tmp_path / "Pfam-A.clans.tsv.gz"
    _write_clans(path, gzipped=True)
    terms = list(PfamLabelProvider().iter_terms(path))
    assert {t.accession for t in terms} == {"PF00001", "PF00008", "PF00069"}


def test_load_writes_labels_and_normalizes_version() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        # Versioned meta should match unversioned dump via normalize.
        _seed_meta(
            conn,
            ["PF00001.15", "PF00069", "not-a-pfam", "PF999999"],
        )
        conn.close()

        clans = work / "fixture.tsv"
        _write_clans(clans)
        stats = load_namespace_labels(
            work, "pfam", source_path=clans, force_download=False
        )
        # PF00001.15 → PF00001, PF00069; not-a-pfam stays dirty; PF999999 missing
        assert stats.labels_written == 2
        assert stats.meta_missing_label == 2  # not-a-pfam + PF999999
        assert stats.fts_written == stats.labels_written + stats.aliases_written

        conn = sqlite3.connect(str(db))
        labels = {
            (acc, lab)
            for acc, lab in conn.execute(
                "SELECT accession, label FROM xref_label WHERE namespace='pfam'"
            )
        }
        assert (
            "PF00001",
            "7 transmembrane receptor (rhodopsin family)",
        ) in labels
        assert ("PF00069", "Protein kinase domain") in labels
        aliases = {
            alias
            for (alias,) in conn.execute(
                "SELECT alias FROM xref_label_alias "
                "WHERE namespace='pfam' AND accession='PF00001'"
            )
        }
        assert aliases == {"7tm_1"}
        assert "GPCR_A" not in aliases
        conn.close()
