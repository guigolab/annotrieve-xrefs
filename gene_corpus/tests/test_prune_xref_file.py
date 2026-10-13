"""Tests for prune_xref_file and searchable-only index writes."""
from __future__ import annotations

import sqlite3
from pathlib import Path

from gene_corpus.namespaces.tiers import SEARCHABLE_NAMESPACES
from gene_corpus.sql.per_gff import (
    connect_per_gff,
    create_indexes,
    init_schema,
    insert_genes,
    insert_xrefs,
    shard_sqlite_path,
    write_shard_meta,
)
from gene_corpus.sql.prune_xref_file import apply_prune, summarize_prune
from gene_corpus.sql.schema import connect_for_build, init_schema as init_global
from gene_corpus.sql.xref_index import (
    ensure_schema,
    index_annotation,
)
from gene_corpus.namespaces.tiers import TIER_A


def _seed_global(db_path: Path) -> sqlite3.Connection:
    conn = connect_for_build(db_path)
    init_global(conn)
    ensure_schema(conn)
    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (1, 'ann-mixed', 9606, 'GCA_1', 'human', 'GenBank',
                  NULL, 'genbank', NULL),
                 (2, 'ann-demoted-only', 9606, 'GCA_2', 'human', 'Ensembl',
                  NULL, 'ensembl', NULL)
        """
    )
    conn.commit()
    return conn


def test_summarize_and_apply_prune(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    conn = _seed_global(db)
    conn.executemany(
        "INSERT INTO xref_file VALUES (?, ?, ?)",
        [
            ("symbol", "tp53", 1),
            ("ensembl", "ENSG1", 1),
            ("ensembl", "ENSG2", 2),
            ("uniprot/trembl", "P1", 2),
        ],
    )
    conn.executemany(
        "INSERT INTO xref_file_status VALUES (?, ?)",
        [(1, 1), (2, 1)],
    )
    conn.commit()

    before = summarize_prune(conn)
    assert before.rows_total == 4
    assert before.rows_to_delete == 3
    assert before.rows_kept == 1
    assert before.placed_to_zero == 1  # key 2 has only demoted
    assert ("ensembl", 2) in before.namespaces_to_delete

    # Dry-run path does not mutate
    assert (
        conn.execute("SELECT COUNT(*) FROM xref_file").fetchone()[0] == 4
    )

    result = apply_prune(conn, db_path=db)
    assert result.rows_deleted == 3
    assert result.placed_flipped == 1
    assert result.vacuumed is True

    rows = conn.execute(
        "SELECT namespace, accession, annotation_key FROM xref_file "
        "ORDER BY namespace, accession"
    ).fetchall()
    assert rows == [("symbol", "tp53", 1)]
    assert (
        conn.execute(
            "SELECT placed FROM xref_file_status WHERE annotation_key = 1"
        ).fetchone()[0]
        == 1
    )
    assert (
        conn.execute(
            "SELECT placed FROM xref_file_status WHERE annotation_key = 2"
        ).fetchone()[0]
        == 0
    )

    after = summarize_prune(conn)
    assert after.rows_to_delete == 0
    assert after.placed_to_zero == 0
    conn.close()


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
        [(1, ns, acc, None, "attr", 0) for ns, acc in xrefs],
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


def test_index_annotation_skips_demoted_namespaces(tmp_path: Path) -> None:
    per_gff = tmp_path / "per_gff"
    db = tmp_path / "gene_corpus.sqlite"
    conn = connect_for_build(db)
    init_global(conn)
    ensure_schema(conn)
    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (1, 'ann-mixed', 9606, 'GCA_1', 'human', 'GenBank',
                  NULL, 'genbank', NULL),
                 (2, 'ann-ensembl-only', 9606, 'GCA_2', 'human', 'Ensembl',
                  NULL, 'ensembl', NULL)
        """
    )
    conn.commit()

    _write_shard(
        per_gff,
        "ann-mixed",
        xrefs=[("symbol", "tp53"), ("ensembl", "ENSG1")],
    )
    _write_shard(
        per_gff,
        "ann-ensembl-only",
        xrefs=[("ensembl", "ENSG9")],
    )

    assert index_annotation(conn, per_gff, 1, "ann-mixed") is True
    rows = {
        (r[0], r[1])
        for r in conn.execute(
            "SELECT namespace, accession FROM xref_file "
            "WHERE annotation_key = 1"
        )
    }
    assert rows == {("symbol", "tp53")}
    assert "ensembl" not in SEARCHABLE_NAMESPACES

    assert index_annotation(conn, per_gff, 2, "ann-ensembl-only") is False
    assert (
        conn.execute(
            "SELECT placed FROM xref_file_status WHERE annotation_key = 2"
        ).fetchone()[0]
        == 0
    )
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM xref_file WHERE annotation_key = 2"
        ).fetchone()[0]
        == 0
    )
    conn.close()
