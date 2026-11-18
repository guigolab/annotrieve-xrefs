"""Tests for the sidecar-based loci build."""
from __future__ import annotations

import sqlite3

import pytest

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.sql import loci_from_sidecars as lfs
from gene_corpus.sql.sidecar import sidecar_path_for
from gene_corpus.tests.test_presence_from_sidecars import corpus


def published(tmp_path, n=5):
    work = corpus(tmp_path, n=n)
    assert serve_build_main(["build", "--work-dir", str(work)]) == 0
    return work


def loci_db(work):
    return sqlite3.connect(work / "serve" / "current" / "loci.sqlite")


def test_loci_rows_match_genes_and_presence(tmp_path):
    work = published(tmp_path)
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 0
    pres = sqlite3.connect(work / "serve" / "current" / "presence.sqlite")
    names = dict(pres.execute("SELECT name, name_id FROM name_dict"))
    ords = dict(pres.execute("SELECT annotation_id, ord FROM annotation"))
    c = loci_db(work)

    meta = dict(c.execute("SELECT key, value FROM meta"))
    assert meta["data_version"] == dict(pres.execute("SELECT key, value FROM meta"))["data_version"]

    rows = c.execute("SELECT ord, local_id, start, end, strand, via, gene_id, primary_name "
                     "FROM locus WHERE name_id=? AND local_id=1 ORDER BY ord", (names["tp53"],)).fetchall()
    assert [r[0] for r in rows] == sorted(ords[f"ann{i}"] for i in range(5))   # tp53 symbol in all 5
    assert all(r[1] == 1 and (r[2], r[3]) == (100, 150) and r[5] == 0 for r in rows)
    assert rows[0][6] == "gene-1"
    # alias-only gene: ann odd -> gene 2 is alias of tp53 too, but gene 2 is not the symbol gene
    alias = c.execute("SELECT ord, local_id, via FROM locus WHERE name_id=? AND local_id=2",
                      (names["tp53"],)).fetchall()
    assert sorted(a[0] for a in alias) == sorted(ords[f"ann{i}"] for i in range(5) if i % 2)
    assert all(a[2] == 1 for a in alias)
    assert c.execute("SELECT n_loci FROM name_loci WHERE name_id=?", (names["tp53"],)).fetchone()[0] \
        == 5 + 2
    n_rows = c.execute("SELECT COUNT(*) FROM locus").fetchone()[0]
    assert n_rows == int(meta["n_loci"]) == sum(r[0] for r in c.execute("SELECT n_loci FROM name_loci"))
    assert set(c.execute("SELECT name FROM dict_feature_type")) == {("gene",)}
    assert dict(c.execute("SELECT seqid_id, seqid FROM seqid WHERE ord=?", (ords["ann0"],))) == {1: "chr1"}


@pytest.mark.parametrize("budget,limit", [(1, 150_000), (300_000, 1)])
def test_batching_and_scratch_store_give_identical_rows(tmp_path, monkeypatch, budget, limit):
    work = published(tmp_path)
    serve_build_main(["loci", "--work-dir", str(work)])
    ref = sorted(loci_db(work).execute("SELECT * FROM locus"))
    monkeypatch.setattr(lfs, "DICT_GENES_LIMIT", limit)          # limit=1 forces the SQLite store
    assert serve_build_main(["loci", "--work-dir", str(work), "--gene-budget", str(budget)]) == 0
    assert sorted(loci_db(work).execute("SELECT * FROM locus")) == ref


def test_refuses_partial_presence_and_changed_manifest(tmp_path, capsys):
    work = corpus(tmp_path)
    serve_build_main(["build", "--work-dir", str(work), "--name-prefix", "uniq"])   # partial only
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 2                 # nothing published
    serve_build_main(["build", "--work-dir", str(work)])
    m = sqlite3.connect(work / "serve" / "manifest.sqlite")
    m.execute("UPDATE annotation SET in_catalogue=0 WHERE annotation_id='ann0'")
    m.commit()
    capsys.readouterr()
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 1
    assert "annotation set differs" in capsys.readouterr().err
    assert not (work / "serve" / "current" / "loci.sqlite").exists()


def test_corrupt_sidecar_fails_without_publishing(tmp_path):
    work = published(tmp_path)
    side = sidecar_path_for(work / "per_gff" / "ann2" / "genes.sqlite")
    side.write_bytes(side.read_bytes()[:-20])
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 1
    assert not (work / "serve" / "current" / "loci.sqlite").exists()
