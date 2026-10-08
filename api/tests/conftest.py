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


def _init_global(path: Path) -> None:
    conn = sqlite3.connect(str(path))
    conn.executescript(
        """
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT UNIQUE NOT NULL,
            taxid INTEGER NOT NULL,
            assembly_accession TEXT NOT NULL,
            organism_name TEXT,
            source_database TEXT NOT NULL,
            source_provider TEXT,
            profile_id TEXT NOT NULL,
            gff_path TEXT
        );
        CREATE TABLE annotation_lineage (
            annotation_key INTEGER NOT NULL,
            taxid INTEGER NOT NULL,
            PRIMARY KEY (annotation_key, taxid)
        ) WITHOUT ROWID;
        CREATE INDEX annotation_lineage_taxid
            ON annotation_lineage (taxid, annotation_key);
        CREATE TABLE gene_hit (
            namespace TEXT NOT NULL,
            accession TEXT NOT NULL,
            annotation_key INTEGER NOT NULL,
            local_id INTEGER NOT NULL,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand INTEGER NOT NULL,
            feature_type TEXT NOT NULL,
            biotype TEXT,
            primary_name TEXT,
            PRIMARY KEY (namespace, accession, annotation_key, local_id)
        ) WITHOUT ROWID;
        CREATE TABLE xref_meta (
            namespace TEXT NOT NULL,
            accession TEXT NOT NULL,
            n_annotations INTEGER NOT NULL,
            n_loci INTEGER NOT NULL,
            PRIMARY KEY (namespace, accession)
        ) WITHOUT ROWID;
        CREATE INDEX xref_meta_ns_loci
            ON xref_meta (namespace, n_loci DESC, accession);
        CREATE TABLE namespace_stats (
            namespace TEXT PRIMARY KEY,
            accession_count INTEGER NOT NULL
        ) WITHOUT ROWID;
        """
    )
    # Human annotation (taxid 9606, lineage includes Mammalia 40674)
    conn.execute(
        """
        INSERT INTO annotation VALUES
        (1, 'ann-human', 9606, 'GCF_000001405.40', 'Homo sapiens',
         'RefSeq', 'NCBI', 'refseq', '/tmp/h.gff.gz')
        """
    )
    # Mouse annotation
    conn.execute(
        """
        INSERT INTO annotation VALUES
        (2, 'ann-mouse', 10090, 'GCF_000001635.27', 'Mus musculus',
         'RefSeq', 'NCBI', 'refseq', '/tmp/m.gff.gz')
        """
    )
    for ak, taxids in (
        (1, (9606, 9605, 40674, 33208)),
        (2, (10090, 10088, 40674, 33208)),
    ):
        for t in taxids:
            conn.execute(
                "INSERT INTO annotation_lineage VALUES (?, ?)", (ak, t)
            )

    # symbol:tp53 on human; alias:tp53 on human; symbol:tp53 on mouse
    hits = [
        ("symbol", "tp53", 1, 10, "NC_000017.11", 7661779, 7687550, -1,
         "gene", "protein_coding", "TP53"),
        ("alias", "tp53", 1, 10, "NC_000017.11", 7661779, 7687550, -1,
         "gene", "protein_coding", "TP53"),
        ("symbol", "tp53", 2, 3, "NC_000011.7", 69580359, 69590348, 1,
         "gene", "protein_coding", "Trp53"),
        ("go", "GO:0008150", 1, 10, "NC_000017.11", 7661779, 7687550, -1,
         "gene", "protein_coding", "TP53"),
        ("go", "GO:0008150", 2, 3, "NC_000011.7", 69580359, 69590348, 1,
         "gene", "protein_coding", "Trp53"),
        ("interpro", "IPR002117", 1, 10, "NC_000017.11", 7661779, 7687550, -1,
         "gene", "protein_coding", "TP53"),
    ]
    conn.executemany(
        """
        INSERT INTO gene_hit VALUES
        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        hits,
    )
    conn.executemany(
        "INSERT INTO xref_meta VALUES (?, ?, ?, ?)",
        [
            ("symbol", "tp53", 2, 2),
            ("alias", "tp53", 1, 1),
            ("go", "GO:0008150", 2, 2),
            ("interpro", "IPR002117", 1, 1),
        ],
    )
    conn.executemany(
        "INSERT INTO namespace_stats VALUES (?, ?)",
        [
            ("alias", 1),
            ("go", 1),
            ("interpro", 1),
            ("symbol", 1),
        ],
    )
    conn.commit()
    conn.close()


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
            (20, "symbol", "brca1", "BRCA1", "attr", 0),
        ],
    )
    conn.commit()
    conn.close()


@pytest.fixture()
def corpus_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "gene_corpus"
    root.mkdir()
    _init_global(root / "gene_corpus.sqlite")
    _init_shard(root / "per_gff" / "ann-human" / "genes.sqlite")
    monkeypatch.setenv("GENE_CORPUS_DIR", str(root))
    # Reload settings cache used by gene_corpus_dir.
    import settings as settings_mod
    import helpers.session as session_mod

    settings_mod.settings.GENE_CORPUS_DIR = str(root)
    session_mod.close_reader()
    yield root
    session_mod.close_reader()


@pytest.fixture()
def client(corpus_dir: Path) -> TestClient:
    from main import create_app

    app = create_app()
    with TestClient(app) as test_client:
        yield test_client
