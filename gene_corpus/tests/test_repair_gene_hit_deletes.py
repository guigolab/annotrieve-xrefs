"""Tests for gene_hit symbol/alias repair via xref_file anti-join."""
from __future__ import annotations

from pathlib import Path

import pytest

from gene_corpus.repair_gene_hit_deletes_cli import main as cli_main
from gene_corpus.sql.repair_gene_hit_deletes import (
    XrefFileGuardError,
    check_xref_file_guard,
    cleanup_spill,
    finalize,
    is_finalized,
    is_scan_done,
    part_tsv_path,
    read_finalize_cursor,
    run_scan,
    scan_complete,
    write_finalize_cursor,
    write_scan_done,
)
from gene_corpus.sql.schema import connect_for_build, init_schema as init_global
from gene_corpus.sql.xref_index import ensure_schema


def _hit(
    ns: str,
    acc: str,
    key: int,
    local_id: int,
    *,
    name: str = "GPX4",
) -> tuple:
    return (
        ns,
        acc,
        key,
        local_id,
        "chr1",
        local_id * 100,
        local_id * 100 + 50,
        1,
        "gene",
        None,
        name,
    )


def _seed(db: Path) -> None:
    """
    Seed a small corpus (stale fraction kept well under 50%):

    - ann 1 indexed: valid symbols gpx4/brca1/tp53 (+extra loci), stale gpx4-201
    - ann 2 indexed: valid aliases phgpx/a1/a2/a3, stale junkalias
    - ann 3 NOT in xref_file_status: gene_hit symbol orphan (must not delete)
    - geneid / go always kept
    """
    conn = connect_for_build(db)
    init_global(conn)
    ensure_schema(conn)
    conn.executemany(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (?, ?, 9606, 'GCA_x', 'human', 'Ensembl', '', 'ensembl', NULL)
        """,
        [(1, "ann-a"), (2, "ann-b"), (3, "ann-unindexed")],
    )
    hits = [
        # valid symbols on ann-a (enough that 2 stale / total < 0.5)
        _hit("symbol", "gpx4", 1, 1),
        _hit("symbol", "brca1", 1, 1),
        _hit("symbol", "tp53", 1, 1),
        _hit("symbol", "egfr", 1, 1),
        _hit("symbol", "myc", 1, 1),
        # stale symbol on ann-a (missing from xref_file)
        _hit("symbol", "gpx4-201", 1, 1),
        _hit("symbol", "gpx4-201", 1, 2),
        # valid aliases on ann-b
        _hit("alias", "phgpx", 2, 1),
        _hit("alias", "a1", 2, 1),
        _hit("alias", "a2", 2, 1),
        _hit("alias", "a3", 2, 1),
        # stale alias on ann-b
        _hit("alias", "junkalias", 2, 1),
        # other namespace kept even if not in xref_file
        _hit("geneid", "12345", 1, 1),
        _hit("go", "GO:0000001", 1, 1),
        # unindexed annotation: stale-looking symbol must NOT be deleted
        _hit("symbol", "unindexed-sym", 3, 1, name="X"),
    ]
    conn.executemany(
        """
        INSERT INTO gene_hit (
            namespace, accession, annotation_key, local_id,
            seqid, start, end, strand, feature_type, biotype, primary_name
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        hits,
    )
    conn.executemany(
        "INSERT INTO xref_file VALUES (?, ?, ?)",
        [
            ("symbol", "gpx4", 1),
            ("symbol", "brca1", 1),
            ("symbol", "tp53", 1),
            ("symbol", "egfr", 1),
            ("symbol", "myc", 1),
            ("alias", "phgpx", 2),
            ("alias", "a1", 2),
            ("alias", "a2", 2),
            ("alias", "a3", 2),
            ("go", "GO:0000001", 1),
            ("geneid", "12345", 1),
        ],
    )
    conn.executemany(
        "INSERT INTO xref_file_status VALUES (?, ?)",
        [(1, 1), (2, 1)],
    )
    conn.commit()
    conn.close()


def _seed_high_stale(db: Path) -> None:
    """Most symbol rows stale → fraction guard should trip."""
    conn = connect_for_build(db)
    init_global(conn)
    ensure_schema(conn)
    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (1, 'ann-a', 9606, 'GCA_x', 'human', 'Ensembl',
                  '', 'ensembl', NULL)
        """
    )
    conn.executemany(
        """
        INSERT INTO gene_hit (
            namespace, accession, annotation_key, local_id,
            seqid, start, end, strand, feature_type, biotype, primary_name
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            _hit("symbol", "keep", 1, 1),
            _hit("symbol", "stale1", 1, 1),
            _hit("symbol", "stale2", 1, 1),
            _hit("alias", "keep-a", 1, 1),
        ],
    )
    conn.executemany(
        "INSERT INTO xref_file VALUES (?, ?, ?)",
        [
            ("symbol", "keep", 1),
            ("alias", "keep-a", 1),
        ],
    )
    conn.execute("INSERT INTO xref_file_status VALUES (1, 1)")
    conn.commit()
    conn.close()


def test_full_apply_deletes_stale_only(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)

    summary = run_scan(
        conn, db_path=db, spill_dir=spill, apply=True, force=False
    )
    assert scan_complete(spill)
    assert summary.stale_rows == 3  # 2x gpx4-201 + junkalias
    by_ns = {n.namespace: n for n in summary.namespaces}
    assert by_ns["symbol"].stale_rows == 2
    assert by_ns["alias"].stale_rows == 1
    assert "gpx4-201" in by_ns["symbol"].sample_accessions
    assert "junkalias" in by_ns["alias"].sample_accessions

    text_s = part_tsv_path(spill, "symbol").read_text(encoding="utf-8")
    text_a = part_tsv_path(spill, "alias").read_text(encoding="utf-8")
    assert "gpx4-201" in text_s
    assert "junkalias" in text_a
    assert "gpx4\t" not in text_s.replace("gpx4-201", "")
    assert "phgpx" not in text_a
    assert "unindexed-sym" not in text_s
    assert "geneid" not in text_s and "geneid" not in text_a

    fin = finalize(conn, spill_dir=spill)
    assert fin.deleted == 3
    assert is_finalized(spill)

    remaining = {
        (r[0], r[1], r[2], r[3])
        for r in conn.execute(
            "SELECT namespace, accession, annotation_key, local_id "
            "FROM gene_hit"
        )
    }
    assert ("symbol", "gpx4", 1, 1) in remaining
    assert ("symbol", "brca1", 1, 1) in remaining
    assert ("alias", "phgpx", 2, 1) in remaining
    assert ("geneid", "12345", 1, 1) in remaining
    assert ("go", "GO:0000001", 1, 1) in remaining
    assert ("symbol", "unindexed-sym", 3, 1) in remaining
    assert ("symbol", "gpx4-201", 1, 1) not in remaining
    assert ("symbol", "gpx4-201", 1, 2) not in remaining
    assert ("alias", "junkalias", 2, 1) not in remaining

    # xref_file untouched
    assert (
        conn.execute("SELECT COUNT(*) FROM xref_file").fetchone()[0] == 11
    )
    cleanup_spill(spill)
    assert not spill.exists()
    conn.close()


def test_dry_run_no_spill_no_deletes(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)
    before = conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0]
    summary = run_scan(
        conn, db_path=db, spill_dir=spill, apply=False, force=False
    )
    assert summary.stale_rows == 3
    assert not spill.exists()
    assert (
        conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0] == before
    )
    conn.close()


def test_status_guard_skips_unindexed(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)
    run_scan(conn, db_path=db, spill_dir=spill, apply=True)
    finalize(conn, spill_dir=spill)
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM gene_hit "
            "WHERE accession = 'unindexed-sym'"
        ).fetchone()[0]
        == 1
    )
    conn.close()


def test_xref_file_guard_aborts_without_force(tmp_path: Path) -> None:
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
        ) VALUES (1, 'ann-a', 9606, 'GCA_x', 'human', 'Ensembl',
                  '', 'ensembl', NULL)
        """
    )
    conn.execute(
        """
        INSERT INTO gene_hit (
            namespace, accession, annotation_key, local_id,
            seqid, start, end, strand, feature_type, biotype, primary_name
        ) VALUES ('symbol', 'gpx4', 1, 1, 'chr1', 1, 100, 1, 'gene', NULL, 'X')
        """
    )
    conn.execute("INSERT INTO xref_file_status VALUES (1, 1)")
    # xref_file empty for symbol
    conn.commit()

    with pytest.raises(XrefFileGuardError, match="xref_file"):
        check_xref_file_guard(conn, force=False)

    # force allows
    counts = check_xref_file_guard(conn, force=True)
    assert counts["symbol"][0] == 1
    assert counts["symbol"][1] == 0

    spill = tmp_path / "spill"
    summary = run_scan(
        conn, db_path=db, spill_dir=spill, apply=True, force=True
    )
    assert summary.stale_rows == 1
    fin = finalize(conn, spill_dir=spill)
    assert fin.deleted == 1
    assert (
        conn.execute("SELECT COUNT(*) FROM gene_hit").fetchone()[0] == 0
    )
    conn.close()


def test_resume_mid_finalize(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)
    run_scan(conn, db_path=db, spill_dir=spill, apply=True)

    # Pretend we already deleted the alias row; leave cursor at that PK.
    conn.execute(
        """
        DELETE FROM gene_hit
        WHERE namespace = 'alias' AND accession = 'junkalias'
          AND annotation_key = 2 AND local_id = 1
        """
    )
    conn.commit()
    write_finalize_cursor(spill, ("alias", "junkalias", 2, 1))

    fin = finalize(conn, spill_dir=spill)
    assert fin.skipped_resume >= 1
    assert fin.deleted == 2  # remaining gpx4-201 rows
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM gene_hit WHERE accession = 'gpx4-201'"
        ).fetchone()[0]
        == 0
    )
    assert is_finalized(spill)
    conn.close()


def test_double_finalize_is_noop(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)
    run_scan(conn, db_path=db, spill_dir=spill, apply=True)
    first = finalize(conn, spill_dir=spill)
    assert first.deleted == 3
    after = set(
        conn.execute(
            "SELECT namespace, accession, annotation_key, local_id "
            "FROM gene_hit"
        ).fetchall()
    )
    # Spill still present for second call (don't cleanup yet)
    second = finalize(conn, spill_dir=spill)
    assert second.deleted == 0
    assert second.tsv_lines == 0
    assert (
        set(
            conn.execute(
                "SELECT namespace, accession, annotation_key, local_id "
                "FROM gene_hit"
            ).fetchall()
        )
        == after
    )
    conn.close()


def test_idempotent_rescan_after_done(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed(db)
    conn = connect_for_build(db)
    first = run_scan(conn, db_path=db, spill_dir=spill, apply=True)
    text0 = part_tsv_path(spill, "symbol").read_text(encoding="utf-8")
    second = run_scan(conn, db_path=db, spill_dir=spill, apply=True)
    text1 = part_tsv_path(spill, "symbol").read_text(encoding="utf-8")
    assert first.stale_rows == second.stale_rows == 3
    assert text0 == text1
    assert scan_complete(spill)
    conn.close()


def test_finalize_requires_scan_complete(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    spill.mkdir()
    _seed(db)
    conn = connect_for_build(db)
    write_scan_done(spill, "alias")
    # symbol done missing
    with pytest.raises(RuntimeError, match="incomplete"):
        finalize(conn, spill_dir=spill)
    conn.close()


def test_read_finalize_cursor_roundtrip(tmp_path: Path) -> None:
    spill = tmp_path / "spill"
    spill.mkdir()
    row = ("symbol", "gpx4-201", 1, 2)
    write_finalize_cursor(spill, row)
    assert read_finalize_cursor(spill) == row


def test_crash_mid_scan_rewrites_partial_tsv(tmp_path: Path) -> None:
    """Partial spill without .done is rewritten on the next apply scan."""
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    spill.mkdir()
    _seed(db)
    # Simulate crash: partial TSV, no done marker.
    part_tsv_path(spill, "symbol").write_text(
        "symbol\tgpx4-201\t1\t1\n", encoding="utf-8"
    )
    assert not is_scan_done(spill, "symbol")

    conn = connect_for_build(db)
    summary = run_scan(conn, db_path=db, spill_dir=spill, apply=True)
    assert scan_complete(spill)
    text = part_tsv_path(spill, "symbol").read_text(encoding="utf-8")
    # Full rewrite includes both gpx4-201 loci, not the truncated crash file.
    assert text.count("gpx4-201") == 2
    assert summary.stale_rows == 3
    conn.close()


def test_stale_fraction_aborts_without_force(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    spill = tmp_path / "spill"
    _seed_high_stale(db)
    conn = connect_for_build(db)

    with pytest.raises(XrefFileGuardError, match="stale fraction"):
        run_scan(conn, db_path=db, spill_dir=spill, apply=False, force=False)

    # No done markers after a rejected apply scan.
    with pytest.raises(XrefFileGuardError, match="stale fraction"):
        run_scan(conn, db_path=db, spill_dir=spill, apply=True, force=False)
    assert not scan_complete(spill)

    summary = run_scan(
        conn, db_path=db, spill_dir=spill, apply=True, force=True
    )
    assert summary.stale_rows == 2
    assert scan_complete(spill)
    fin = finalize(conn, spill_dir=spill)
    assert fin.deleted == 2
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM gene_hit WHERE namespace = 'symbol'"
        ).fetchone()[0]
        == 1
    )
    conn.close()


def test_cli_finalize_only_requires_apply(tmp_path: Path) -> None:
    db = tmp_path / "gene_corpus.sqlite"
    _seed(db)
    rc = cli_main(
        ["--corpus-db", str(db), "--finalize-only", "--spill-dir", str(tmp_path / "spill")]
    )
    assert rc == 1
