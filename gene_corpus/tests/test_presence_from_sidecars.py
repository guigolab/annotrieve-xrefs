"""Tests for the sidecar-based presence build (build plan M4)."""
from __future__ import annotations

import sqlite3

import pytest

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.sql.postings import decode_postings
from gene_corpus.sql.presence_from_sidecars import build_from_sidecars
from gene_corpus.sql.serve_presence import BuildError
from gene_corpus.sql.sidecar import sidecar_path_for
from gene_corpus.tests.test_catalogue_backfill import tsv
from gene_corpus.tests.test_sidecar_manifest import make_shard

def genes(*ids):
    return [(i, f"gene-{i}", "gene", "chr1", i * 100, i * 100 + 50, 1, None, None, None, None, 0)
            for i in ids]


def corpus(tmp_path, n=5):
    work = tmp_path / "corpus"
    root = work / "per_gff"
    for i in range(n):
        xrefs = [(1, "symbol", "tp53", None, "gff", 0), (2, "symbol", f"uniq{i}", None, "gff", 0),
                 (2, "alias", "p53" if i % 2 == 0 else "tp53", None, "gff", 0)]
        make_shard(root, f"ann{i}", genes=genes(1, 2), xrefs=xrefs)
    make_shard(root, "nonames", genes=genes(1), xrefs=[])
    make_shard(root, "stale", genes=genes(1), xrefs=[(1, "symbol", "tp53", None, "gff", 0)])
    ids = [f"ann{i}" for i in range(n)] + ["nonames"]       # "stale" is not in the catalogue
    rep = tmp_path / "report.tsv"
    rep.write_text("\n".join(tsv(*[(a, "GCF_1", "n", "o", 9606, "RefSeq", "p", "2024-01-01", "1;9606")
                                   for a in ids])) + "\n")
    assert serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"]) == 0
    assert serve_build_main(["catalogue", "--work-dir", str(work), "--tsv-file", str(rep)]) == 0
    return work


def build(work, tmp_path, **kw):
    return build_from_sidecars(work / "serve" / "manifest.sqlite", work / "per_gff",
                               tmp_path / "presence.sqlite", tmp_dir=tmp_path / "tmp",
                               log=lambda m: None, **kw)


def postings(path, name):
    c = sqlite3.connect(path)
    sym, alias = c.execute("SELECT p.sym, p.alias FROM name_dict d JOIN posting p USING (name_id) "
                           "WHERE d.name=?", (name,)).fetchone()
    return (decode_postings(sym) if sym else []), (decode_postings(alias) if alias else [])


@pytest.mark.parametrize("batch", [1, 2, 256])
def test_build_postings_and_filters(tmp_path, batch):
    work = corpus(tmp_path)
    stats = build(work, tmp_path, batch=batch)
    out = tmp_path / "presence.sqlite"
    c = sqlite3.connect(out)
    ords = dict(c.execute("SELECT annotation_id, ord FROM annotation"))
    assert set(ords) == {f"ann{i}" for i in range(5)} | {"nonames"}      # stale shard excluded
    assert stats["n_annotations_with_names"] == 5
    sym, alias = postings(out, "tp53")
    assert sym == sorted(ords[f"ann{i}"] for i in range(5))
    assert alias == sorted(ords[f"ann{i}"] for i in range(5) if i % 2)   # odd ones alias tp53 too
    assert postings(out, "uniq3") == ([ords["ann3"]], [])
    assert postings(out, "p53") == ([], sorted(ords[f"ann{i}"] for i in range(5) if i % 2 == 0))
    n_alias_only = c.execute("SELECT n_alias_only FROM name_dict WHERE name='p53'").fetchone()[0]
    assert n_alias_only == 3
    assert c.execute("SELECT COUNT(*) FROM lineage WHERE taxid=9606").fetchone()[0] == 6
    assert stats["peak_rss_mib"] < 500


def test_build_is_deterministic_across_batch_sizes(tmp_path):
    work = corpus(tmp_path)
    v1 = build(work, tmp_path, batch=1)["data_version"]
    v2 = build(work, tmp_path, batch=256)["data_version"]
    assert v1 == v2
    rows = lambda: sorted(sqlite3.connect(tmp_path / "presence.sqlite").execute(  # noqa: E731
        "SELECT d.name, d.n_sym, p.sym, p.alias FROM name_dict d JOIN posting p USING (name_id)"))
    a = rows()
    build(work, tmp_path, batch=2)
    assert rows() == a


def test_stale_sidecar_or_manifest_mismatch_fails(tmp_path):
    work = corpus(tmp_path)
    m = sqlite3.connect(work / "serve" / "manifest.sqlite")
    m.execute("UPDATE annotation SET n_name_rows = n_name_rows + 1 WHERE annotation_id='ann1'")
    m.commit()
    with pytest.raises(BuildError, match="n_name_rows"):
        build(work, tmp_path)
    m.execute("UPDATE annotation SET shard_fingerprint='v1:0:0:0:0' WHERE annotation_id='ann1'")
    m.commit()
    with pytest.raises(BuildError, match="stale"):
        build(work, tmp_path)


def test_truncated_sidecar_fails(tmp_path):
    work = corpus(tmp_path)
    side = sidecar_path_for(work / "per_gff" / "ann2" / "genes.sqlite")
    side.write_bytes(side.read_bytes()[:-20])
    with pytest.raises(Exception):
        build(work, tmp_path)


def test_prefix_trial_and_cli_publish(tmp_path, capsys):
    work = corpus(tmp_path)
    assert serve_build_main(["build", "--work-dir", str(work), "--name-prefix", "uniq"]) == 0
    assert "partial build written" in capsys.readouterr().out
    assert not (work / "serve" / "current").exists()
    assert serve_build_main(["build", "--work-dir", str(work)]) == 0
    assert (work / "serve" / "current" / "presence.sqlite").is_file()
    assert not any((work / "serve" / "tmp").iterdir())          # runs are cleaned
