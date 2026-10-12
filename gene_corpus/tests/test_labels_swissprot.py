"""Swiss-Prot FASTA header parse + load against xref_meta."""
from __future__ import annotations

import gzip
import sqlite3
import tempfile
from pathlib import Path

from gene_corpus.labels.load import load_namespace_labels
from gene_corpus.labels.providers.swissprot import (
    SwissprotLabelProvider,
    parse_swissprot_header,
)
from gene_corpus.sql.schema import connect_for_build, init_schema

_FASTA_FIXTURE = """\
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
>sp|P10978|GAG_DROME Gag polyprotein OS=Drosophila melanogaster OX=7227 GN=gag PE=1 SV=2
MGQKISNITSSSTSNSTSVSTTSTSSTSTSTSTSTSTSTSTST
>sp|P04323|ENV_DROME Envelope polyprotein OS=Drosophila melanogaster OX=7227 PE=1 SV=1
MKSLLLLLLL
"""


def _write_fasta(path: Path, *, gzipped: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = _FASTA_FIXTURE.encode("utf-8")
    if gzipped:
        path.write_bytes(gzip.compress(data))
    else:
        path.write_bytes(data)


def _seed_meta(conn: sqlite3.Connection, accessions: list[str]) -> None:
    conn.executemany(
        "INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci) "
        "VALUES ('uniprot/swiss-prot', ?, 1, 1)",
        [(a,) for a in accessions],
    )
    conn.commit()


def test_parse_header_maps_ac_name_entry_gn() -> None:
    term = parse_swissprot_header(
        ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 "
        "OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4\n"
    )
    assert term is not None
    assert term.accession == "P04637"
    assert term.label == "Cellular tumor antigen p53"
    assert "P53_HUMAN" in term.aliases
    assert "TP53" in term.aliases


def test_parse_header_without_gn() -> None:
    term = parse_swissprot_header(
        ">sp|P04323|ENV_DROME Envelope polyprotein "
        "OS=Drosophila melanogaster OX=7227 PE=1 SV=1\n"
    )
    assert term is not None
    assert term.accession == "P04323"
    assert term.label == "Envelope polyprotein"
    assert term.aliases == ("ENV_DROME",)


def test_iter_terms_gzipped_fasta(tmp_path: Path) -> None:
    path = tmp_path / "uniprot_sprot.fasta.gz"
    _write_fasta(path, gzipped=True)
    terms = list(SwissprotLabelProvider().iter_terms(path))
    by_id = {t.accession: t for t in terms}
    assert set(by_id) == {"P04637", "P10978", "P04323"}
    assert by_id["P10978"].label == "Gag polyprotein"
    assert "gag" in by_id["P10978"].aliases


def test_load_writes_labels_for_meta_ids() -> None:
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        _seed_meta(conn, ["P04637", "P10978", "P99999"])
        conn.close()

        fasta = work / "fixture.fasta"
        _write_fasta(fasta, gzipped=False)
        stats = load_namespace_labels(
            work, "uniprot/swiss-prot", source_path=fasta, force_download=False
        )
        assert stats.labels_written == 2
        assert stats.meta_missing_label == 1

        conn = sqlite3.connect(str(db))
        labels = {
            (acc, lab)
            for acc, lab in conn.execute(
                "SELECT accession, label FROM xref_label "
                "WHERE namespace='uniprot/swiss-prot'"
            )
        }
        assert ("P04637", "Cellular tumor antigen p53") in labels
        assert ("P10978", "Gag polyprotein") in labels
        n_fts = conn.execute(
            "SELECT COUNT(*) FROM xref_label_fts "
            "WHERE namespace='uniprot/swiss-prot'"
        ).fetchone()[0]
        assert n_fts == stats.fts_written
        conn.close()

        # Idempotent reload
        stats2 = load_namespace_labels(
            work, "uniprot/swiss-prot", source_path=fasta, force_download=False
        )
        assert stats2.labels_written == 2
        conn = sqlite3.connect(str(db))
        n2 = conn.execute(
            "SELECT COUNT(*) FROM xref_label "
            "WHERE namespace='uniprot/swiss-prot'"
        ).fetchone()[0]
        assert n2 == 2
        conn.close()
