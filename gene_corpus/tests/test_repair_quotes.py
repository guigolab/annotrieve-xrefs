"""Repair quote-wrapped symbol/alias accessions in gene_corpus.sqlite."""
from __future__ import annotations

import tempfile
from pathlib import Path

from gene_corpus.repair_cli import main as repair_main
from gene_corpus.sql.repair_quotes import (
    AFFECTED_BARE_STATE_NAME,
    count_wrapped_accessions,
    repair_wrapped_accessions,
)
from gene_corpus.sql.schema import connect_for_build, init_schema, seed_annotations
from gene_corpus.stream.symbol_gate import normalize_symbol, strip_wrapping_quotes


def test_strip_wrapping_quotes_and_normalize() -> None:
    assert strip_wrapping_quotes("'TP53'") == "TP53"
    assert strip_wrapping_quotes('"TP53"') == "TP53"
    assert strip_wrapping_quotes("TP53") == "TP53"
    assert strip_wrapping_quotes("'x") == "'x"
    assert normalize_symbol("'TP53'") == "tp53"
    assert normalize_symbol("TP53") == "tp53"
    assert normalize_symbol("  'Abc'  ") == "abc"


def _seed_db(path: Path) -> None:
    conn = connect_for_build(path)
    init_schema(conn)
    seed_annotations(
        conn,
        [
            {
                "annotation_key": 1,
                "annotation_id": "ann1",
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
                "annotation_id": "ann2",
                "taxid": 10090,
                "assembly_accession": "GCA_2",
                "organism_name": "mouse",
                "source_database": "GenBank",
                "source_provider": None,
                "profile_id": "genbank",
                "gff_path": None,
            },
        ],
    )
    # Collision: same locus under quoted + bare
    conn.execute(
        """
        INSERT INTO gene_hit (
            namespace, accession, annotation_key, local_id,
            seqid, start, end, strand, feature_type, biotype, primary_name
        ) VALUES
            ('symbol', '''a1cf''', 1, 10, 'chr1', 1, 10, 1, 'gene', NULL, 'A1CF'),
            ('symbol', 'a1cf', 1, 10, 'chr1', 1, 10, 1, 'gene', NULL, 'A1CF'),
            ('symbol', '''a1cf''', 2, 3, 'chr2', 5, 15, -1, 'gene', NULL, 'A1cf'),
            ('symbol', '''onlyq''', 1, 20, 'chr3', 1, 2, 1, 'gene', NULL, 'ONLYQ')
        """
    )
    conn.execute(
        """
        INSERT INTO xref_meta (namespace, accession, n_annotations, n_loci)
        VALUES
            ('symbol', '''a1cf''', 2, 2),
            ('symbol', 'a1cf', 1, 1),
            ('symbol', '''onlyq''', 1, 1),
            ('alias', 'p53', 1, 1)
        """
    )
    conn.execute(
        """
        INSERT INTO namespace_stats (namespace, accession_count)
        VALUES ('symbol', 3), ('alias', 1)
        """
    )
    conn.commit()
    conn.close()


def test_repair_collision_and_rename() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        db = root / "gene_corpus.sqlite"
        state = root / "repair_quotes"
        _seed_db(db)

        conn = connect_for_build(db)
        result = repair_wrapped_accessions(
            conn, state_dir=state, apply=True
        )
        assert result.applied
        assert result.quoted_keys == 2
        assert result.dup_loci_skipped == 1  # ann1/local 10 already bare
        assert result.moved_rows == 2  # ann2 a1cf + onlyq
        assert result.wrapped_remaining == 0

        rows = {
            (r[0], r[1], r[2], r[3])
            for r in conn.execute(
                "SELECT namespace, accession, annotation_key, local_id "
                "FROM gene_hit ORDER BY 1,2,3,4"
            )
        }
        assert ("symbol", "'a1cf'", 1, 10) not in rows
        assert ("symbol", "'onlyq'", 1, 20) not in rows
        assert ("symbol", "a1cf", 1, 10) in rows
        assert ("symbol", "a1cf", 2, 3) in rows
        assert ("symbol", "onlyq", 1, 20) in rows

        meta = {
            (r[0], r[1]): (r[2], r[3])
            for r in conn.execute(
                "SELECT namespace, accession, n_annotations, n_loci FROM xref_meta"
            )
        }
        assert ("symbol", "'a1cf'") not in meta
        assert ("symbol", "'onlyq'") not in meta
        assert meta[("symbol", "a1cf")] == (2, 2)
        assert meta[("symbol", "onlyq")] == (1, 1)

        stats = dict(
            conn.execute(
                "SELECT namespace, accession_count FROM namespace_stats"
            )
        )
        assert stats["symbol"] == 2  # a1cf + onlyq
        assert stats["alias"] == 1
        conn.close()


def test_repair_dry_run_unchanged() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        db = root / "gene_corpus.sqlite"
        state = root / "repair_quotes"
        _seed_db(db)

        conn = connect_for_build(db)
        before_hit = conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0]
        before_meta = list(
            conn.execute(
                "SELECT namespace, accession, n_annotations, n_loci "
                "FROM xref_meta ORDER BY 1,2"
            )
        )
        result = repair_wrapped_accessions(
            conn, state_dir=state, apply=False
        )
        assert not result.applied
        assert result.quoted_keys == 2
        assert result.moved_rows == 0
        assert conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0] == before_hit
        after_meta = list(
            conn.execute(
                "SELECT namespace, accession, n_annotations, n_loci "
                "FROM xref_meta ORDER BY 1,2"
            )
        )
        assert after_meta == before_meta
        assert (state / "quoted.txt").is_file()
        assert (state / AFFECTED_BARE_STATE_NAME).is_file()
        assert count_wrapped_accessions(conn, ("symbol", "alias")) == 2
        conn.close()


def test_repair_cli_recompute_only() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        db = root / "gene_corpus.sqlite"
        state = root / "repair_quotes"
        _seed_db(db)

        # Simulate interrupted apply: moves done, meta stale for bare a1cf
        conn = connect_for_build(db)
        repair_wrapped_accessions(conn, state_dir=state, apply=True)
        # Corrupt meta on purpose
        conn.execute(
            "UPDATE xref_meta SET n_annotations = 99, n_loci = 99 "
            "WHERE namespace = 'symbol' AND accession = 'a1cf'"
        )
        conn.commit()
        conn.close()

        rc = repair_main(
            ["--work-dir", str(root), "--state-dir", str(state), "--recompute-only"]
        )
        assert rc == 0

        conn = connect_for_build(db)
        n_ann, n_loci = conn.execute(
            "SELECT n_annotations, n_loci FROM xref_meta "
            "WHERE namespace = 'symbol' AND accession = 'a1cf'"
        ).fetchone()
        assert (n_ann, n_loci) == (2, 2)
        conn.close()


def test_repair_cli_apply() -> None:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        _seed_db(root / "gene_corpus.sqlite")
        rc = repair_main(["--work-dir", str(root), "--apply"])
        assert rc == 0
        conn = connect_for_build(root / "gene_corpus.sqlite")
        assert count_wrapped_accessions(conn, ("symbol",)) == 0
        conn.close()
