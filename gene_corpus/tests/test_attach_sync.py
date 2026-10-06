"""Incremental attach + sync lineage refresh tests."""
from __future__ import annotations

import tempfile
from pathlib import Path

from gene_corpus.sql.attach import (
    attach_shard,
    next_annotation_key,
    read_seed_row_from_shard,
    rebuild_namespace_stats,
)
from gene_corpus.sql.lineage import (
    TAXONOMY_SHA256_META_KEY,
    get_meta_value,
    load_parent_map,
    refresh_lineage_if_taxonomy_changed,
    taxonomy_file_sha256,
)
from gene_corpus.sql.per_gff import (
    connect_per_gff,
    init_schema as init_per_gff_schema,
    shard_dir,
    shard_sqlite_path,
    write_shard_meta,
)
from gene_corpus.sql.schema import connect_for_build, init_schema
from gene_corpus.sync_cli import main as sync_main


def _write_taxonomy_tsv(path: Path) -> Path:
    path.write_text(
        "taxid\tparent_taxid\tscientific_name\n"
        "9606\t9605\tHomo sapiens\n"
        "9605\t40674\tHomo\n"
        "40674\t\tMammalia\n"
        "10090\t10088\tMus musculus\n"
        "10088\t40674\tMus\n",
        encoding="utf-8",
    )
    return path


def _seed_merge_shard(
    per_gff: Path,
    annotation_id: str,
    *,
    genes: list[tuple[int, str, int, int]],
    xrefs: list[tuple[int, str, str]],
    tier_a: list[tuple[str, str, int]],
    taxid: str = "9606",
) -> Path:
    folder = shard_dir(per_gff, annotation_id)
    folder.mkdir(parents=True, exist_ok=True)
    path = shard_sqlite_path(per_gff, annotation_id)
    conn = connect_per_gff(path)
    init_per_gff_schema(conn)
    for local_id, seqid, start, end in genes:
        conn.execute(
            """
            INSERT INTO gene (
                local_id, source_gene_id, feature_type, seqid, start, end,
                strand, biotype, primary_name, prose, prose_kind, parse_flags
            ) VALUES (?, ?, 'gene', ?, ?, ?, 1, 'protein_coding', ?, NULL, NULL, 0)
            """,
            (local_id, f"g{local_id}", seqid, start, end, f"gene{local_id}"),
        )
    for local_id, namespace, accession in xrefs:
        conn.execute(
            """
            INSERT INTO gene_xref (
                local_id, namespace, accession, display, origin, via_level
            ) VALUES (?, ?, ?, NULL, 'attr', 0)
            """,
            (local_id, namespace, accession),
        )
    conn.executemany(
        """
        INSERT INTO tier_a_counts (namespace, accession, gene_count)
        VALUES (?, ?, ?)
        """,
        tier_a,
    )
    write_shard_meta(
        conn,
        {
            "annotation_id": annotation_id,
            "taxid": taxid,
            "assembly_accession": "GCA",
            "organism_name": "",
            "source_database": "GenBank",
            "source_provider": "",
            "profile_id": "genbank",
            "gff_path": "",
        },
    )
    conn.commit()
    conn.close()
    return path


def test_attach_assigns_increasing_keys() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        per_gff = root / "per_gff"
        tax_tsv = _write_taxonomy_tsv(root / "taxonomy.tsv")
        parents = load_parent_map(tax_tsv)
        _seed_merge_shard(
            per_gff,
            "a",
            genes=[(0, "chr1", 1, 10)],
            xrefs=[(0, "symbol", "tp53")],
            tier_a=[("symbol", "tp53", 1)],
            taxid="9606",
        )
        _seed_merge_shard(
            per_gff,
            "b",
            genes=[(0, "chr2", 1, 5)],
            xrefs=[(0, "go", "GO:0008150")],
            tier_a=[("go", "GO:0008150", 1)],
            taxid="10090",
        )
        db = root / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)

        key_a = next_annotation_key(conn)
        assert key_a == 1
        seed_a = read_seed_row_from_shard(per_gff, "a", key_a)
        assert attach_shard(conn, per_gff, seed_a, taxonomy_parents=parents)

        key_b = next_annotation_key(conn)
        assert key_b == 2
        seed_b = read_seed_row_from_shard(per_gff, "b", key_b)
        assert attach_shard(conn, per_gff, seed_b, taxonomy_parents=parents)

        keys = conn.execute(
            "SELECT annotation_key, annotation_id FROM annotation "
            "ORDER BY annotation_key"
        ).fetchall()
        assert keys == [(1, "a"), (2, "b")]

        # Idempotent re-attach
        assert not attach_shard(
            conn, per_gff, seed_a, taxonomy_parents=parents
        )
        assert next_annotation_key(conn) == 3

        rebuild_namespace_stats(conn)
        meta = conn.execute(
            "SELECT namespace, accession, n_annotations, n_loci FROM xref_meta "
            "ORDER BY namespace, accession"
        ).fetchall()
        assert meta == [
            ("go", "GO:0008150", 1, 1),
            ("symbol", "tp53", 1, 1),
        ]
        stats = conn.execute(
            "SELECT namespace, accession_count FROM namespace_stats "
            "ORDER BY namespace"
        ).fetchall()
        assert stats == [("go", 1), ("symbol", 1)]
        hits = conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0]
        assert hits == 2
        lineage_a = conn.execute(
            "SELECT taxid FROM annotation_lineage "
            "WHERE annotation_key = 1 ORDER BY taxid"
        ).fetchall()
        assert lineage_a == [(9605,), (9606,), (40674,)]
        conn.close()


def test_refresh_lineage_on_taxonomy_hash_change() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        per_gff = root / "per_gff"
        tax_tsv = _write_taxonomy_tsv(root / "taxonomy.tsv")
        parents = load_parent_map(tax_tsv)
        _seed_merge_shard(
            per_gff,
            "a",
            genes=[(0, "chr1", 1, 10)],
            xrefs=[(0, "symbol", "tp53")],
            tier_a=[("symbol", "tp53", 1)],
        )
        db = root / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        seed = read_seed_row_from_shard(per_gff, "a", 1)
        assert attach_shard(conn, per_gff, seed, taxonomy_parents=parents)

        refreshed, n1 = refresh_lineage_if_taxonomy_changed(
            conn, tax_tsv, parents
        )
        assert refreshed is True
        assert n1 >= 3
        digest = taxonomy_file_sha256(tax_tsv)
        assert get_meta_value(conn, TAXONOMY_SHA256_META_KEY) == digest

        refreshed2, n2 = refresh_lineage_if_taxonomy_changed(
            conn, tax_tsv, parents
        )
        assert refreshed2 is False
        assert n2 == n1

        # Parent change for 9606 → rebuild
        tax_tsv.write_text(
            "taxid\tparent_taxid\tscientific_name\n"
            "9606\t9999\tHomo sapiens\n"
            "9999\t\tNewParent\n",
            encoding="utf-8",
        )
        parents2 = load_parent_map(tax_tsv)
        refreshed3, n3 = refresh_lineage_if_taxonomy_changed(
            conn, tax_tsv, parents2
        )
        assert refreshed3 is True
        lineage = conn.execute(
            "SELECT taxid FROM annotation_lineage "
            "WHERE annotation_key = 1 ORDER BY taxid"
        ).fetchall()
        assert lineage == [(9606,), (9999,)]
        assert n3 == 2
        conn.close()


def test_attach_failure_rolls_back_annotation() -> None:
    """Corrupt shard must not leave a dangling annotation row."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        per_gff = root / "per_gff"
        tax_tsv = _write_taxonomy_tsv(root / "taxonomy.tsv")
        parents = load_parent_map(tax_tsv)
        folder = shard_dir(per_gff, "bad")
        folder.mkdir(parents=True)
        shard = shard_sqlite_path(per_gff, "bad")
        shard.write_bytes(b"not a sqlite database")

        db = root / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        seed = {
            "annotation_key": 1,
            "annotation_id": "bad",
            "taxid": 9606,
            "assembly_accession": "GCA",
            "organism_name": None,
            "source_database": "GenBank",
            "source_provider": None,
            "profile_id": "genbank",
            "gff_path": None,
        }
        try:
            attach_shard(conn, per_gff, seed, taxonomy_parents=parents)
            raise AssertionError("expected attach to fail")
        except Exception:
            pass
        n = conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
        assert n == 0
        assert next_annotation_key(conn) == 1
        conn.close()


def test_seed_rejects_meta_id_mismatch() -> None:
    with tempfile.TemporaryDirectory() as d:
        per_gff = Path(d) / "per_gff"
        _seed_merge_shard(
            per_gff,
            "dir-id",
            genes=[(0, "chr1", 1, 10)],
            xrefs=[(0, "symbol", "tp53")],
            tier_a=[("symbol", "tp53", 1)],
        )
        # Overwrite meta annotation_id so it disagrees with the directory name.
        shard = shard_sqlite_path(per_gff, "dir-id")
        conn = connect_per_gff(shard)
        write_shard_meta(
            conn,
            {
                "annotation_id": "other-id",
                "taxid": "9606",
                "assembly_accession": "GCA",
                "organism_name": "",
                "source_database": "GenBank",
                "source_provider": "",
                "profile_id": "genbank",
                "gff_path": "",
            },
        )
        conn.commit()
        conn.close()
        try:
            read_seed_row_from_shard(per_gff, "dir-id", 1)
            raise AssertionError("expected ValueError")
        except ValueError as exc:
            assert "does not match" in str(exc)


def test_sync_cli_attaches_new_shard_from_local_report() -> None:
    """Sync with local --report / --taxonomy-tsv (no network)."""
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        per_gff = work / "per_gff"
        files_root = work / "files"
        files_root.mkdir()
        tax_tsv = _write_taxonomy_tsv(work / "taxonomy.tsv")

        # Bootstrap DB empty of annotations
        db = work / "gene_corpus.sqlite"
        conn = connect_for_build(db)
        init_schema(conn)
        conn.close()

        # Pre-built shard for annotation "new1" (skip real GFF harvest)
        _seed_merge_shard(
            per_gff,
            "new1",
            genes=[(0, "chr1", 1, 10)],
            xrefs=[(0, "symbol", "brca1")],
            tier_a=[("symbol", "brca1", 1)],
            taxid="9606",
        )

        report = work / "report.tsv"
        report.write_text(
            "annotation_id\tassembly_accession\torganism_name\ttaxid\t"
            "database\tprovider\tsource_url\tbgzip_path\t"
            "attribute_keys\tfeature_types\thas_cds\troot_type_counts\n"
            "new1\tGCA_TEST\tHomo sapiens\t9606\t"
            "GenBank\t\t\tmissing.gff.gz\t"
            "ID;Name;Dbxref\tgene;CDS\ttrue\t{\"gene\":1}\n",
            encoding="utf-8",
        )

        rc = sync_main(
            [
                "--work-dir",
                str(work),
                "--files-root",
                str(files_root),
                "--report",
                str(report),
                "--taxonomy-tsv",
                str(tax_tsv),
            ]
        )
        # resolve skips missing GFF; shard already on disk → attach still runs
        assert rc == 0
        conn = connect_for_build(db)
        keys = conn.execute(
            "SELECT annotation_key, annotation_id FROM annotation"
        ).fetchall()
        assert keys == [(1, "new1")]
        n_hit = conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0]
        assert n_hit == 1
        assert get_meta_value(conn, TAXONOMY_SHA256_META_KEY) is not None
        conn.close()
