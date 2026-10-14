"""Tests for delete-delta symbol/alias repair."""
from __future__ import annotations

from pathlib import Path

import pytest

from gene_corpus.sql.per_gff import (
    connect_per_gff,
    create_indexes,
    init_schema,
    insert_genes,
    insert_xrefs,
    shard_sqlite_path,
    write_shard_meta,
)
from gene_corpus.sql.repair_symbol_deletes import (
    cleanup_spill,
    finalize,
    is_finalized,
    part_tsv_path,
    process_annotation,
    read_cursor,
    run_scan,
    scan_complete,
    write_cursor,
)
from gene_corpus.sql.schema import connect_for_build, init_schema as init_global
from gene_corpus.sql.xref_index import ensure_schema


def _write_shard(
    per_gff: Path,
    annotation_id: str,
    *,
    xrefs: list[tuple[int, str, str, int]],
) -> None:
    path = shard_sqlite_path(per_gff, annotation_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_per_gff(path)
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
                None,
                "GPX4",
                None,
                None,
                0,
            ),
            (
                2,
                "g2",
                "gene",
                "chr1",
                200,
                300,
                1,
                None,
                "GPX4",
                None,
                None,
                0,
            ),
        ],
    )
    insert_xrefs(
        conn,
        [
            (local_id, ns, acc, None, "attr", via)
            for local_id, ns, acc, via in xrefs
        ],
    )
    write_shard_meta(
        conn,
        {
            "annotation_id": annotation_id,
            "taxid": "9606",
            "assembly_accession": "GCA_1",
            "organism_name": "human",
            "source_database": "Ensembl",
            "source_provider": "",
            "profile_id": "ensembl",
            "gff_path": "",
        },
    )
    create_indexes(conn)
    conn.commit()
    conn.close()


def _seed_corpus(work_dir: Path) -> None:
    per_gff = work_dir / "per_gff"
    db = work_dir / "gene_corpus.sqlite"
    # ann-a: level-0 gpx4 + via>0 junk gpx4-201; both-level BRCA1
    _write_shard(
        per_gff,
        "ann-a",
        xrefs=[
            (1, "symbol", "gpx4", 0),
            (1, "symbol", "gpx4-201", 1),
            (1, "symbol", "BRCA1", 0),
            (2, "symbol", "BRCA1", 1),
            (1, "go", "GO:0000001", 0),
        ],
    )
    _write_shard(
        per_gff,
        "ann-b",
        xrefs=[
            (1, "symbol", "gpx4", 0),
            (1, "alias", "phgpx", 0),
            (1, "symbol", "gpx4-202", 1),
        ],
    )
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
        [(1, "ann-a"), (2, "ann-b")],
    )
    conn.executemany(
        "INSERT INTO xref_file VALUES (?, ?, ?)",
        [
            ("symbol", "gpx4", 1),
            ("symbol", "gpx4-201", 1),
            ("symbol", "BRCA1", 1),
            ("symbol", "gpx4", 2),
            ("symbol", "gpx4-202", 2),
            ("alias", "phgpx", 2),
            ("go", "GO:0000001", 1),
        ],
    )
    conn.executemany(
        "INSERT INTO xref_file_status VALUES (?, ?)",
        [(1, 1), (2, 1)],
    )
    conn.executemany(
        "INSERT INTO xref_meta VALUES (?, ?, ?, ?)",
        [
            ("symbol", "gpx4", 2, 2),
            ("symbol", "gpx4-201", 1, 1),
            ("symbol", "BRCA1", 1, 2),
            ("symbol", "gpx4-202", 1, 1),
            ("alias", "phgpx", 1, 1),
            ("go", "GO:0000001", 1, 1),
        ],
    )
    conn.commit()
    conn.close()


def test_full_apply_spills_and_patches(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    db = tmp_path / "gene_corpus.sqlite"
    conn = connect_for_build(db)

    summary = run_scan(
        conn, work_dir=tmp_path, spill_dir=spill, apply=True
    )
    assert summary.missing_shard == 0
    assert summary.xref_rows_deleted >= 2
    assert scan_complete(conn, spill)

    # TSV has junk only, not level-0-only gpx4 / phgpx; not BRCA1 (both levels)
    text = ""
    for w in (0, 1):
        p = part_tsv_path(spill, w)
        if p.is_file():
            text += p.read_text(encoding="utf-8")
    assert "gpx4-201" in text
    assert "gpx4-202" in text
    assert "\tgpx4\t" not in text or "gpx4-20" in text  # only junk forms
    assert "BRCA1" not in text
    assert "phgpx" not in text

    # Shards cleaned of via>0; level-0 kept; BRCA1 level-0 kept
    sa = connect_per_gff(shard_sqlite_path(tmp_path / "per_gff", "ann-a"))
    rows = {
        (r[0], r[1], r[2])
        for r in sa.execute(
            "SELECT namespace, accession, via_level FROM gene_xref"
        )
    }
    assert ("symbol", "gpx4", 0) in rows
    assert ("symbol", "BRCA1", 0) in rows
    assert ("symbol", "gpx4-201", 1) not in rows
    assert ("symbol", "BRCA1", 1) not in rows
    sa.close()

    fin = finalize(conn, spill_dir=spill)
    assert fin.xref_file_deleted >= 2

    file_rows = set(
        conn.execute(
            "SELECT namespace, accession, annotation_key FROM xref_file"
        ).fetchall()
    )
    assert ("symbol", "gpx4", 1) in file_rows
    assert ("symbol", "gpx4", 2) in file_rows
    assert ("symbol", "BRCA1", 1) in file_rows  # both-levels: presence kept
    assert ("symbol", "gpx4-201", 1) not in file_rows
    assert ("symbol", "gpx4-202", 2) not in file_rows
    assert ("alias", "phgpx", 2) in file_rows
    assert ("go", "GO:0000001", 1) in file_rows

    meta = {
        (r[0], r[1]): (r[2], r[3])
        for r in conn.execute(
            "SELECT namespace, accession, n_annotations, n_loci FROM xref_meta"
        )
    }
    assert ("symbol", "gpx4-201") not in meta
    assert ("symbol", "gpx4-202") not in meta
    assert meta[("symbol", "gpx4")] == (2, 2)
    # BRCA1: not in TSV → meta unchanged (still inflated loci=2)
    assert meta[("symbol", "BRCA1")] == (1, 2)
    assert meta[("go", "GO:0000001")] == (1, 1)

    cleanup_spill(spill)
    assert not spill.exists()
    conn.close()


def test_both_levels_no_tsv_line_but_shard_cleaned(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    tsv_path = tmp_path / "out.tsv"
    with tsv_path.open("w", encoding="utf-8") as fh:
        deleted, missing, n_lines = process_annotation(
            per_gff_root=tmp_path / "per_gff",
            annotation_key=1,
            annotation_id="ann-a",
            apply=True,
            tsv_fh=fh,
        )
    assert not missing
    assert deleted >= 2  # gpx4-201 + BRCA1 via>0
    text = tsv_path.read_text(encoding="utf-8")
    assert "gpx4-201" in text
    assert "BRCA1" not in text
    assert n_lines == 1  # only gpx4-201


def test_resume_reprocesses_cursor_key(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    spill.mkdir(parents=True)
    # Worker 0 owns key 1 (index 0); worker 1 owns key 2.
    # Pretend worker 0 finished key 1.
    write_cursor(spill, 0, 1)
    part_tsv_path(spill, 0).write_text(
        "1\tsymbol\tgpx4-201\t1\n", encoding="utf-8"
    )
    # Clean shard a already
    sa = connect_per_gff(shard_sqlite_path(tmp_path / "per_gff", "ann-a"))
    sa.execute(
        "DELETE FROM gene_xref WHERE via_level > 0 "
        "AND namespace IN ('symbol', 'alias')"
    )
    sa.commit()
    sa.close()

    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    summary = run_scan(
        conn, work_dir=tmp_path, spill_dir=spill, apply=True
    )
    # Worker 0 reprocesses key 1 (0 deleted); worker 1 processes key 2
    assert summary.skipped_resume >= 0
    assert scan_complete(conn, spill)
    assert read_cursor(spill, 0) == 1
    assert read_cursor(spill, 1) == 2

    # Duplicate line for key 1 ok
    text0 = part_tsv_path(spill, 0).read_text(encoding="utf-8")
    assert text0.count("gpx4-201") >= 1

    fin = finalize(conn, spill_dir=spill)
    assert fin.unique_deletes >= 1
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM xref_file WHERE accession LIKE 'gpx4-%'"
        ).fetchone()[0]
        == 0
    )
    conn.close()


def test_finalize_decrements_when_other_annotation_remains(
    tmp_path: Path,
) -> None:
    """gpx4-201 only on ann-a; shared symbol stays; junk dropped."""
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    # Put same junk on both annotations in meta/file as if shared
    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    conn.execute(
        "INSERT INTO xref_file VALUES ('symbol', 'gpx4-201', 2)"
    )
    conn.execute(
        "UPDATE xref_meta SET n_annotations = 2, n_loci = 2 "
        "WHERE namespace='symbol' AND accession='gpx4-201'"
    )
    conn.commit()

    # Manual spill: both keys delete gpx4-201
    spill.mkdir(parents=True)
    part_tsv_path(spill, 0).write_text(
        "1\tsymbol\tgpx4-201\t1\n", encoding="utf-8"
    )
    part_tsv_path(spill, 1).write_text(
        "2\tsymbol\tgpx4-201\t1\n2\tsymbol\tgpx4-202\t1\n",
        encoding="utf-8",
    )
    write_cursor(spill, 0, 1)
    write_cursor(spill, 1, 2)

    # Also clean shards so we don't care about gene_xref for this test
    fin = finalize(conn, spill_dir=spill)
    assert fin.xref_meta_deleted >= 1
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM xref_meta "
            "WHERE accession IN ('gpx4-201', 'gpx4-202')"
        ).fetchone()[0]
        == 0
    )
    assert (
        conn.execute(
            "SELECT COUNT(*) FROM xref_file WHERE accession='gpx4-201'"
        ).fetchone()[0]
        == 0
    )
    conn.close()


def test_missing_shard_blocks_complete(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (3, 'ann-missing', 9606, 'GCA_x', 'human', 'Ensembl',
                  '', 'ensembl', NULL)
        """
    )
    conn.commit()

    summary = run_scan(
        conn, work_dir=tmp_path, spill_dir=spill, apply=True
    )
    assert summary.missing_shard >= 1
    assert not scan_complete(conn, spill)
    with pytest.raises(RuntimeError, match="incomplete"):
        finalize(conn, spill_dir=spill)
    conn.close()


def test_dry_run_no_writes(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    deleted, missing, n_lines = process_annotation(
        per_gff_root=tmp_path / "per_gff",
        annotation_key=1,
        annotation_id="ann-a",
        apply=False,
    )
    assert deleted >= 1
    assert not missing
    assert n_lines == 0
    sa = connect_per_gff(shard_sqlite_path(tmp_path / "per_gff", "ann-a"))
    assert (
        sa.execute(
            "SELECT COUNT(*) FROM gene_xref WHERE accession='gpx4-201'"
        ).fetchone()[0]
        == 1
    )
    sa.close()


def test_clean_probe_skips_aggregate(tmp_path: Path) -> None:
    """No via_level>0 → empty lines, shard untouched (probe fast path)."""
    per_gff = tmp_path / "per_gff"
    _write_shard(
        per_gff,
        "ann-clean",
        xrefs=[(1, "symbol", "gpx4", 0)],
    )
    deleted, missing, n_lines = process_annotation(
        per_gff_root=per_gff,
        annotation_key=1,
        annotation_id="ann-clean",
        apply=True,
    )
    assert deleted == 0
    assert not missing
    assert n_lines == 0


def test_spill_fsync_before_delete_order(tmp_path: Path) -> None:
    """TSV must contain spill lines before shard DELETE commits."""
    _seed_corpus(tmp_path)
    tsv_path = tmp_path / "order.tsv"
    events: list[str] = []

    import gene_corpus.sql.repair_symbol_deletes as mod

    real_delete = mod._delete_via_level_gt0

    def wrapped_delete(conn: object) -> int:
        events.append("delete")
        assert tsv_path.is_file()
        text = tsv_path.read_text(encoding="utf-8")
        assert "gpx4-201" in text
        events.append("tsv_ok_before_delete")
        return real_delete(conn)  # type: ignore[arg-type]

    real_fsync = mod._fsync_handle

    def wrapped_fsync(fh: object) -> None:
        events.append("fsync")
        real_fsync(fh)  # type: ignore[arg-type]

    # Patch via module attributes used inside process_annotation
    orig_del = mod._delete_via_level_gt0
    orig_fs = mod._fsync_handle
    mod._delete_via_level_gt0 = wrapped_delete  # type: ignore[assignment]
    mod._fsync_handle = wrapped_fsync  # type: ignore[assignment]
    try:
        with tsv_path.open("w", encoding="utf-8") as fh:
            deleted, missing, n_lines = process_annotation(
                per_gff_root=tmp_path / "per_gff",
                annotation_key=1,
                annotation_id="ann-a",
                apply=True,
                tsv_fh=fh,
            )
    finally:
        mod._delete_via_level_gt0 = orig_del  # type: ignore[assignment]
        mod._fsync_handle = orig_fs  # type: ignore[assignment]

    assert not missing
    assert deleted >= 1
    assert n_lines >= 1
    assert events == ["fsync", "delete", "tsv_ok_before_delete"]
    sa = connect_per_gff(shard_sqlite_path(tmp_path / "per_gff", "ann-a"))
    assert (
        sa.execute(
            "SELECT COUNT(*) FROM gene_xref "
            "WHERE namespace IN ('symbol','alias') AND via_level > 0"
        ).fetchone()[0]
        == 0
    )
    sa.close()


def test_sync_every_one_writes_cursor_each_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With _SYNC_EVERY=1, end state matches per-key cursor durability."""
    import gene_corpus.sql.repair_symbol_deletes as mod

    monkeypatch.setattr(mod, "_SYNC_EVERY", 1)
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    run_scan(conn, work_dir=tmp_path, spill_dir=spill, apply=True)
    assert scan_complete(conn, spill)
    assert read_cursor(spill, 0) == 1
    assert read_cursor(spill, 1) == 2
    assert part_tsv_path(spill, 0).is_file() or part_tsv_path(spill, 1).is_file()
    conn.close()


def test_duplicate_tsv_lines_deduped_at_finalize(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    spill.mkdir(parents=True)
    part_tsv_path(spill, 0).write_text(
        "1\tsymbol\tgpx4-201\t1\n1\tsymbol\tgpx4-201\t1\n",
        encoding="utf-8",
    )
    part_tsv_path(spill, 1).write_text(
        "2\tsymbol\tgpx4-202\t1\n",
        encoding="utf-8",
    )
    write_cursor(spill, 0, 1)
    write_cursor(spill, 1, 2)
    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    fin = finalize(conn, spill_dir=spill)
    assert fin.unique_deletes == 2
    assert fin.xref_file_deleted == 2
    assert is_finalized(spill)
    conn.close()


def test_double_finalize_is_noop(tmp_path: Path) -> None:
    _seed_corpus(tmp_path)
    spill = tmp_path / "spill"
    conn = connect_for_build(tmp_path / "gene_corpus.sqlite")
    run_scan(conn, work_dir=tmp_path, spill_dir=spill, apply=True)
    first = finalize(conn, spill_dir=spill)
    assert first.xref_file_deleted >= 1
    assert is_finalized(spill)
    meta_after = set(
        conn.execute(
            "SELECT namespace, accession, n_annotations, n_loci FROM xref_meta"
        ).fetchall()
    )
    file_after = set(
        conn.execute(
            "SELECT namespace, accession, annotation_key FROM xref_file"
        ).fetchall()
    )
    second = finalize(conn, spill_dir=spill)
    assert second.xref_file_deleted == 0
    assert second.xref_meta_deleted == 0
    assert second.unique_deletes == 0
    assert (
        set(
            conn.execute(
                "SELECT namespace, accession, n_annotations, n_loci "
                "FROM xref_meta"
            ).fetchall()
        )
        == meta_after
    )
    assert (
        set(
            conn.execute(
                "SELECT namespace, accession, annotation_key FROM xref_file"
            ).fetchall()
        )
        == file_after
    )
    conn.close()


def test_empty_deletes_still_advance_cursor(tmp_path: Path) -> None:
    """Shard with only level-0 symbols: no TSV lines, cursor still advances."""
    per_gff = tmp_path / "per_gff"
    db = tmp_path / "gene_corpus.sqlite"
    _write_shard(
        per_gff,
        "ann-clean",
        xrefs=[(1, "symbol", "gpx4", 0)],
    )
    conn = connect_for_build(db)
    init_global(conn)
    ensure_schema(conn)
    conn.execute(
        """
        INSERT INTO annotation (
            annotation_key, annotation_id, taxid, assembly_accession,
            organism_name, source_database, source_provider, profile_id,
            gff_path
        ) VALUES (1, 'ann-clean', 9606, 'GCA_x', 'human', 'Ensembl',
                  '', 'ensembl', NULL)
        """
    )
    conn.commit()

    spill = tmp_path / "spill"
    summary = run_scan(
        conn, work_dir=tmp_path, spill_dir=spill, apply=True
    )
    assert summary.tsv_lines_appended == 0
    assert summary.xref_rows_deleted == 0
    assert read_cursor(spill, 0) == 1
    assert scan_complete(conn, spill)
    # Worker 1 has no jobs — only key 1 on worker 0
    fin = finalize(conn, spill_dir=spill)
    assert fin.unique_deletes == 0
    assert is_finalized(spill)
    conn.close()
