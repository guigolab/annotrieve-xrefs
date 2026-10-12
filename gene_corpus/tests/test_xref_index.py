"""Tests for gene_corpus.sql.xref_index / index CLI."""
from __future__ import annotations

import sqlite3
from pathlib import Path

from gene_corpus.sql.per_gff import (
    connect_per_gff,
    create_indexes,
    init_schema,
    insert_genes,
    insert_xrefs,
    shard_sqlite_path,
    write_shard_meta,
)
from gene_corpus.sql.schema import connect_for_build, init_schema as init_global
from gene_corpus.sql.xref_index import (
    ensure_schema,
    index_annotation,
    index_missing_annotations,
)
from gene_corpus.namespaces.tiers import TIER_A


def _seed_global(db_path: Path, rows: list[dict]) -> sqlite3.Connection:
    conn = connect_for_build(db_path)
    init_global(conn)
    ensure_schema(conn)
    conn.executemany(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id, gff_path
        ) VALUES (
            :annotation_key, :annotation_id, :taxid, :assembly_accession,
            :organism_name, :source_database, :source_provider, :profile_id,
            :gff_path
        )
        """,
        rows,
    )
    conn.commit()
    return conn


def _write_shard(
    per_gff: Path,
    annotation_id: str,
    *,
    xrefs: list[tuple[str, str]],
) -> None:
    path = shard_sqlite_path(per_gff, annotation_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".sqlite.tmp")
    conn = connect_per_gff(tmp)
    init_schema(conn)
    insert_genes(
        conn,
        [
            (
                1,
                "g1",
                "gene",
                "chr1",
                1,
                100,
                1,
                "protein_coding",
                "GENE",
                None,
                None,
                0,
            )
        ],
    )
    insert_xrefs(
        conn,
        [
            (1, ns, acc, None, "attr", 0)
            for ns, acc in xrefs
        ],
    )
    from gene_corpus.sql.per_gff import build_tier_a_counts

    build_tier_a_counts(conn, TIER_A)
    write_shard_meta(
        conn,
        {
            "annotation_id": annotation_id,
            "taxid": "9606",
            "assembly_accession": "GCA_1",
            "organism_name": "human",
            "source_database": "GenBank",
            "source_provider": "",
            "profile_id": "genbank",
            "gff_path": "/tmp/x.gff.gz",
        },
    )
    create_indexes(conn)
    conn.close()
    tmp.replace(path)


def test_index_annotation_empty_and_placed(tmp_path: Path) -> None:
    per_gff = tmp_path / "per_gff"
    db = tmp_path / "gene_corpus.sqlite"
    conn = _seed_global(
        db,
        [
            {
                "annotation_key": 1,
                "annotation_id": "ann-empty",
                "taxid": 9606,
                "assembly_accession": "GCA_1",
                "organism_name": "human",
                "source_database": "GenBank",
                "source_provider": None,
                "profile_id": "genbank",
                "gff_path": None,
            },
            {
                "annotation_key": 2,
                "annotation_id": "ann-ok",
                "taxid": 9606,
                "assembly_accession": "GCA_2",
                "organism_name": "human",
                "source_database": "Ensembl",
                "source_provider": None,
                "profile_id": "ensembl",
                "gff_path": None,
            },
        ],
    )
    _write_shard(per_gff, "ann-empty", xrefs=[])
    # Empty xrefs still need a finished shard with empty tier_a_counts.
    # _write_shard with empty xrefs builds empty tier_a_counts.
    _write_shard(per_gff, "ann-ok", xrefs=[("symbol", "tp53")])

    assert index_annotation(conn, per_gff, 1, "ann-empty") is False
    status = conn.execute(
        "SELECT placed FROM xref_file_status WHERE annotation_key = 1"
    ).fetchone()
    assert status is not None and int(status[0]) == 0
    assert (
        conn.execute("SELECT COUNT(*) FROM xref_file").fetchone()[0] == 0
    )

    assert index_annotation(conn, per_gff, 2, "ann-ok") is True
    assert (
        conn.execute(
            "SELECT placed FROM xref_file_status WHERE annotation_key = 2"
        ).fetchone()[0]
        == 1
    )
    rows = conn.execute(
        "SELECT namespace, accession, annotation_key FROM xref_file"
    ).fetchall()
    assert rows == [("symbol", "tp53", 2)]

    # Idempotent: second call does not insert again.
    assert index_annotation(conn, per_gff, 2, "ann-ok") is False
    assert (
        conn.execute("SELECT COUNT(*) FROM xref_file").fetchone()[0] == 1
    )
    conn.close()


def test_index_missing_annotations(tmp_path: Path) -> None:
    per_gff = tmp_path / "per_gff"
    db = tmp_path / "gene_corpus.sqlite"
    conn = _seed_global(
        db,
        [
            {
                "annotation_key": 1,
                "annotation_id": "ann-a",
                "taxid": 1,
                "assembly_accession": "GCA_1",
                "organism_name": None,
                "source_database": "RefSeq",
                "source_provider": None,
                "profile_id": "refseq",
                "gff_path": None,
            },
            {
                "annotation_key": 2,
                "annotation_id": "ann-b",
                "taxid": 1,
                "assembly_accession": "GCA_2",
                "organism_name": None,
                "source_database": "RefSeq",
                "source_provider": None,
                "profile_id": "refseq",
                "gff_path": None,
            },
        ],
    )
    _write_shard(per_gff, "ann-a", xrefs=[("go", "GO:0008150")])
    _write_shard(per_gff, "ann-b", xrefs=[])

    placed, empty, errors = index_missing_annotations(conn, per_gff)
    assert errors == 0
    assert placed == 1
    assert empty == 1
    placed2, empty2, errors2 = index_missing_annotations(conn, per_gff)
    assert (placed2, empty2, errors2) == (0, 0, 0)
    conn.close()
