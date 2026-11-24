"""Fixtures: locked serve schema + one per-GFF shard."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

# Ensure api/ is on sys.path when pytest is run from the repo root.
import sys

API_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = API_ROOT.parent
for path in (str(API_ROOT), str(REPO_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)


def _init_shard(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.executescript(
        """
        CREATE TABLE gene (
            local_id INTEGER PRIMARY KEY,
            source_gene_id TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand INTEGER NOT NULL,
            biotype TEXT,
            primary_name TEXT,
            prose TEXT,
            prose_kind TEXT,
            parse_flags INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE gene_xref (
            local_id INTEGER NOT NULL,
            namespace TEXT NOT NULL,
            accession TEXT NOT NULL,
            display TEXT,
            origin TEXT NOT NULL,
            via_level INTEGER NOT NULL,
            UNIQUE (local_id, namespace, accession)
        );
        CREATE TABLE tier_a_counts (
            namespace TEXT NOT NULL,
            accession TEXT NOT NULL,
            gene_count INTEGER NOT NULL,
            PRIMARY KEY (namespace, accession)
        ) WITHOUT ROWID;
        """
    )
    conn.executemany(
        """
        INSERT INTO gene VALUES
        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                10,
                "gene-TP53",
                "gene",
                "NC_000017.11",
                7661779,
                7687550,
                -1,
                "protein_coding",
                "TP53",
                "tumor protein p53",
                "description",
                0,
            ),
            (
                20,
                "gene-BRCA1",
                "gene",
                "NC_000017.11",
                43044295,
                43125483,
                -1,
                "protein_coding",
                "BRCA1",
                "BRCA1 DNA repair associated",
                "description",
                0,
            ),
        ],
    )
    conn.executemany(
        "INSERT INTO gene_xref VALUES (?, ?, ?, ?, ?, ?)",
        [
            (10, "symbol", "tp53", "TP53", "attr", 0),
            (10, "alias", "tp53", "p53", "attr", 0),
            (10, "go", "GO:0008150", None, "attr", 2),
            (10, "ensembl_transcript", "ENST00000269305", None, "attr", 1),
            (
                10,
                "refseq_transcript",
                "XM_015139305.2",
                "XM_015139305.2",
                "attr",
                1,
            ),
            (
                10,
                "uniprot/swiss-prot",
                "P04637",
                "P53_HUMAN",
                "dbxref",
                0,
            ),
            (20, "symbol", "brca1", "BRCA1", "attr", 0),
        ],
    )
    conn.executemany(
        "INSERT INTO tier_a_counts VALUES (?, ?, ?)",
        [
            ("symbol", "tp53", 1),
            ("symbol", "brca1", 1),
            ("alias", "tp53", 1),
            ("go", "GO:0008150", 1),
        ],
    )
    conn.commit()
    conn.close()


@pytest.fixture()
def corpus_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "gene_corpus"
    root.mkdir()
    _init_shard(root / "per_gff" / "ann-human" / "genes.sqlite")
    monkeypatch.setenv("GENE_CORPUS_DIR", str(root))
    import settings as settings_mod

    settings_mod.settings.GENE_CORPUS_DIR = str(root)
    yield root


@pytest.fixture()
def client(corpus_dir: Path) -> TestClient:
    from main import create_app

    app = create_app()
    with TestClient(app) as test_client:
        yield test_client
