"""Tests for gene_corpus.sql.catalogue and the catalogue/sidecars/status commands (build plan M3)."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.sql.catalogue import (
    CatalogueError,
    apply_catalogue,
    load_catalogue,
    report_url,
)
from gene_corpus.sql.manifest import index_shard, open_manifest
from gene_corpus.sql.sidecar import sidecar_path_for
from gene_corpus.tests.test_sidecar_manifest import make_shard

HEADER = ("annotation_id\tassembly_accession\tassembly_name\torganism_name\ttaxid\tdatabase\t"
          "provider\tsource_url\tbgzip_path\tcsi_path\trelease_date\ttaxon_lineage")


def tsv(*rows: tuple) -> list[str]:
    lines = [HEADER]
    for r in rows:
        aid, acc, an, org, taxid, db, prov, rel, lin = r
        lines.append("\t".join([aid, acc, an, org, str(taxid), db, prov, "http://x", "b/p", "b/p.csi", rel, lin]))
    return lines


def test_report_url_has_the_extended_fields():
    assert report_url("https://h/api/v0/") == (
        "https://h/api/v0/annotations/report?selected_fields=release_date,taxon_lineage")


def test_load_fills_metadata_flags_and_lineage(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "a1")
    make_shard(root, "a2")
    make_shard(root, "only_shard")
    conn = open_manifest(tmp_path / "m.sqlite")
    for d in ("a1", "a2", "only_shard"):
        index_shard(conn, root, d)

    stats = load_catalogue(conn, tsv(
        ("a1", "GCF_1", "GRCh38", "Homo sapiens", 9606, "RefSeq", "NCBI", "2024-04-01", "1;131567;9606;40674"),
        ("a2", "GCF_2", "Mmul", "Macaca", 9544, "Ensembl", "ensembl", "2023-01-01", "9544;9443"),
        ("gone", "GCF_3", "X", "Y", 1, "RefSeq", "NCBI", "2020-01-01", "1"),
    ))
    assert stats.rows == 3 and stats.matched_annotations == 2
    assert stats.catalogue_only == 1 and stats.shard_only == ["only_shard"]
    assert stats.database_mismatch == 1   # a2 shard says RefSeq (helper default), catalogue says Ensembl

    row = conn.execute("SELECT assembly_name, release_date, in_catalogue FROM annotation "
                       "WHERE annotation_id='a1'").fetchone()
    assert row == ("GRCh38", "2024-04-01", 1)
    assert conn.execute("SELECT in_catalogue FROM annotation WHERE annotation_id='only_shard'").fetchone()[0] == 0

    ord_a1 = conn.execute("SELECT ord FROM annotation WHERE annotation_id='a1'").fetchone()[0]
    taxids = {t for (t,) in conn.execute("SELECT taxid FROM lineage WHERE ord=?", (ord_a1,))}
    assert taxids == {1, 131567, 9606, 40674}          # the shard's own taxid 9606 is included once


def test_lineage_includes_own_taxid_even_if_report_omits_it(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "a1")                              # shard taxid 9606
    conn = open_manifest(tmp_path / "m.sqlite")
    index_shard(conn, root, "a1")
    load_catalogue(conn, tsv(("a1", "GCF_1", "n", "o", 9606, "RefSeq", "p", "2024-01-01", "1;2759")))
    assert {t for (t,) in conn.execute("SELECT taxid FROM lineage")} == {1, 2759, 9606}


def test_missing_columns_and_empty_report_are_errors(tmp_path):
    conn = open_manifest(tmp_path / "m.sqlite")
    with pytest.raises(CatalogueError):
        load_catalogue(conn, ["annotation_id\ttaxid", "a\t1"])
    with pytest.raises(CatalogueError):
        load_catalogue(conn, [])
    with pytest.raises(CatalogueError):
        load_catalogue(conn, [HEADER])


def test_failed_load_keeps_previous_catalogue(tmp_path):
    conn = open_manifest(tmp_path / "m.sqlite")
    load_catalogue(conn, tsv(("a1", "GCF_1", "n", "o", 1, "RefSeq", "p", "2024-01-01", "1")))
    with pytest.raises(CatalogueError):
        load_catalogue(conn, [HEADER])
    conn.rollback()
    assert conn.execute("SELECT COUNT(*) FROM catalogue").fetchone()[0] == 1


def test_duplicates_and_blank_ids_are_counted(tmp_path):
    conn = open_manifest(tmp_path / "m.sqlite")
    stats = load_catalogue(conn, tsv(
        ("a1", "GCF_1", "n", "o", 1, "RefSeq", "p", "2024-01-01", "1"),
        ("a1", "GCF_1", "n", "o", 1, "RefSeq", "p", "2024-01-01", "1"),
        ("", "GCF_1", "n", "o", 1, "RefSeq", "p", "2024-01-01", "1"),
    ))
    assert (stats.rows, stats.duplicate_ids, stats.skipped_blank_id) == (3, 1, 1)
    assert conn.execute("SELECT COUNT(*) FROM catalogue").fetchone()[0] == 1


def test_apply_catalogue_picks_up_shards_registered_later(tmp_path):
    root = tmp_path / "per_gff"
    conn = open_manifest(tmp_path / "m.sqlite")
    load_catalogue(conn, tsv(("a1", "GCF_1", "GRCh38", "o", 9606, "RefSeq", "p", "2024-04-01", "1;9606")))
    make_shard(root, "a1")
    index_shard(conn, root, "a1")
    assert conn.execute("SELECT in_catalogue FROM annotation").fetchone()[0] == 0
    apply_catalogue(conn)
    conn.commit()
    assert conn.execute("SELECT assembly_name, in_catalogue FROM annotation").fetchone() == ("GRCh38", 1)


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------
def make_corpus(tmp_path: Path, n: int = 6) -> Path:
    work = tmp_path / "corpus"
    root = work / "per_gff"
    sources = ["RefSeq", "Ensembl", "GenBank", "CommunityRegistry", "RefSeq", "Mystery"]
    for i in range(n):
        make_shard(root, f"ann{i:02d}", source=sources[i % len(sources)])
    return work


def test_sidecars_command_parallel_and_resumable(tmp_path, capsys):
    work = make_corpus(tmp_path)
    assert serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "2"]) == 1  # Mystery fails
    out = capsys.readouterr().out
    assert "FAILED ann05" in out and "Mystery" in out

    conn = open_manifest(work / "serve" / "manifest.sqlite")
    status = dict(conn.execute("SELECT annotation_id, status FROM annotation"))
    assert status == {"ann00": "indexed", "ann01": "indexed", "ann02": "indexed",
                      "ann03": "skipped_community", "ann04": "indexed", "ann05": "failed"}
    assert [r[0] for r in conn.execute("SELECT annotation_id FROM annotation ORDER BY ord")] == sorted(status)
    for d in ("ann00", "ann01", "ann02", "ann04"):
        assert sidecar_path_for(work / "per_gff" / d / "genes.sqlite").is_file()
    for d in ("ann03", "ann05"):
        assert not sidecar_path_for(work / "per_gff" / d / "genes.sqlite").exists()

    side = sidecar_path_for(work / "per_gff" / "ann00" / "genes.sqlite")
    before = side.stat().st_mtime_ns
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"])
    assert side.stat().st_mtime_ns == before               # fresh sidecars are not rewritten
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1", "--force"])
    assert side.stat().st_mtime_ns > before


def test_sidecars_ordinals_are_deterministic_with_many_workers(tmp_path):
    work = make_corpus(tmp_path, n=6)
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "3"])
    conn = open_manifest(work / "serve" / "manifest.sqlite")
    assert dict(conn.execute("SELECT annotation_id, ord FROM annotation")) == {
        f"ann{i:02d}": i + 1 for i in range(6)}


def test_sidecars_sample_and_ids_options(tmp_path):
    work = make_corpus(tmp_path)
    ids = tmp_path / "ids.txt"
    ids.write_text("ann00\nann02\n")
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1", "--ids", str(ids)])
    conn = open_manifest(work / "serve" / "manifest.sqlite")
    assert {r[0] for r in conn.execute("SELECT annotation_id FROM annotation")} == {"ann00", "ann02"}


def test_catalogue_command_from_file_then_status(tmp_path, capsys):
    work = make_corpus(tmp_path)
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"])
    f = tmp_path / "report.tsv"
    f.write_text("\n".join(tsv(
        ("ann00", "GCF_1", "GRCh38", "Homo sapiens", 9606, "RefSeq", "p", "2024-04-01", "1;9606"),
        ("ann01", "GCF_2", "X", "Y", 9544, "Ensembl", "p", "2023-01-01", "9544"),
    )) + "\n")
    capsys.readouterr()
    assert serve_build_main(["catalogue", "--work-dir", str(work), "--tsv-file", str(f)]) == 0
    out = capsys.readouterr().out
    assert "catalogue: 2 rows" in out and "matched=2" in out and "catalogue_only=0" in out

    assert serve_build_main(["status", "--work-dir", str(work)]) == 0
    out = capsys.readouterr().out
    assert "indexed" in out and "skipped_community" in out and "failed" in out
    assert "in_catalogue=2" in out


def test_sidecars_reapplies_stored_catalogue(tmp_path):
    work = make_corpus(tmp_path, n=2)
    f = tmp_path / "report.tsv"
    f.write_text("\n".join(tsv(
        ("ann00", "GCF_1", "GRCh38", "o", 9606, "RefSeq", "p", "2024-04-01", "1;9606"))) + "\n")
    serve_build_main(["catalogue", "--work-dir", str(work), "--tsv-file", str(f)])   # manifest empty yet
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"])
    conn = open_manifest(work / "serve" / "manifest.sqlite")
    assert conn.execute("SELECT assembly_name FROM annotation WHERE annotation_id='ann00'").fetchone()[0] == "GRCh38"


def test_legacy_build_arguments_still_dispatch_to_the_old_build(tmp_path, capsys):
    rc = serve_build_main(["--work-dir", str(tmp_path)])     # no gene_corpus.sqlite there
    assert rc == 2
    assert "gene_corpus.sqlite missing" in capsys.readouterr().err


def _try_alloc(_):
    try:
        bytearray(400 * 1024 * 1024)
        return "allocated"
    except MemoryError:
        return "memory_error"


def test_worker_memory_cap_turns_overuse_into_an_exception():
    import sys
    from concurrent.futures import ProcessPoolExecutor

    from gene_corpus.serve_backfill_cli import _init_worker

    if not sys.platform.startswith("linux"):
        pytest.skip("RLIMIT_AS cap is applied on Linux")
    with ProcessPoolExecutor(max_workers=1, initializer=_init_worker, initargs=(256,)) as ex:
        assert ex.submit(_try_alloc, 0).result() == "memory_error"
    with ProcessPoolExecutor(max_workers=1, initializer=_init_worker, initargs=(0,)) as ex:
        assert ex.submit(_try_alloc, 0).result() == "allocated"


def test_missing_lists_only_harvestable_catalogue_annotations(tmp_path, capsys):
    work = make_corpus(tmp_path, n=2)             # ann00 RefSeq, ann01 Ensembl have shards
    f = tmp_path / "report.tsv"
    f.write_text("\n".join(tsv(
        ("ann00", "GCF_1", "n", "o", 9606, "RefSeq", "p", "2024-01-01", "1"),
        ("new_gb", "GCA_2", "n", "o", 1, "GenBank", "p", "2024-01-01", "1"),
        ("new_ens", "GCA_3", "n", "o", 1, "Ensembl", "p", "2024-01-01", "1"),
        ("new_comm", "GCA_4", "n", "o", 1, "CommunityRegistry", "p", "2024-01-01", "1"),
        ("new_odd", "GCA_5", "n", "o", 1, "Mystery", "p", "2024-01-01", "1"),
    )) + "\n")
    serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"])
    serve_build_main(["catalogue", "--work-dir", str(work), "--tsv-file", str(f)])
    capsys.readouterr()
    assert serve_build_main(["missing", "--work-dir", str(work)]) == 0
    out = capsys.readouterr().out
    assert "without a shard: 4" in out and "harvestable (Ensembl/GenBank/RefSeq): 2" in out
    assert "community (skipped): 1" in out and "unknown source (skipped): 1" in out
    ids = (work / "serve" / "missing_ids.txt").read_text().split()
    assert ids == ["new_ens", "new_gb"]
    serve_build_main(["missing", "--work-dir", str(work), "--limit", "1"])
    assert (work / "serve" / "missing_ids.txt").read_text().split() == ["new_ens"]
