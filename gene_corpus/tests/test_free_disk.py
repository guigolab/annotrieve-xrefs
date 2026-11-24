"""tools/free_disk.py: dry run lists, --apply removes only what production does not use."""
from __future__ import annotations

import importlib.util
import json
import os
import sqlite3
import sys
import time
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "free_disk", Path(__file__).resolve().parents[2] / "tools" / "free_disk.py")
free_disk = importlib.util.module_from_spec(SPEC)
sys.modules["free_disk"] = free_disk
SPEC.loader.exec_module(free_disk)


def _file(path: Path, n: int = 100, age_s: float = 0) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * n)
    if age_s:
        t = time.time() - age_s
        os.utime(path, (t, t))
    return path


def _age(path: Path, age_s: float) -> None:
    t = time.time() - age_s
    os.utime(path, (t, t))


@pytest.fixture()
def tree(tmp_path: Path):
    work, results = tmp_path / "corpus", tmp_path / "results"
    serve = work / "serve"
    for i, v in enumerate(("v1", "v2", "v3")):
        _file(serve / v / "presence.sqlite", 1000)
        _age(serve / v, 1000 * (3 - i))                       # v3 newest
    (serve / "current").symlink_to("v3")
    _file(serve / "partial-abc" / "presence.sqlite", 500)
    _file(serve / ".presence.building.sqlite", 300)
    _file(work / "gene_corpus.sqlite", 4000)
    _file(work / "per_gff" / "ann1" / "genes.sqlite", 200)
    _file(work / "per_gff" / "ann1" / "genes.sqlite.tmp123", 50)
    _file(work / "per_gff" / "gone" / "genes.sqlite", 200)
    _file(work / "per_gff" / "removed" / "genes.sqlite", 200)
    m = sqlite3.connect(serve / "manifest.sqlite")
    m.execute("CREATE TABLE annotation (annotation_id TEXT, status TEXT)")
    m.executemany("INSERT INTO annotation VALUES (?,?)", [("ann1", "indexed"), ("removed", "removed")])
    m.commit()
    m.close()
    now = time.time()
    for name, state, age in (("a" * 20, "done", 3 * 3600), ("b" * 20, "done", 60),
                             ("c" * 20, "failed", 2 * 3600), ("d" * 20, "queued", 9 * 3600),
                             ("e" * 20, "running", 30)):
        d = results / name
        _file(d / "result.fa.gz", 100)
        (d / "status.json").write_text(json.dumps({"state": state, "updated": now - age}))
    _file(results / ("f" * 20) / "result.fa", 100)
    _age(results / ("f" * 20), 2 * 3600)                    # no status file, old
    _file(results / ".quota" / "x.json", 10)
    return work, results


def names(items):
    return sorted((i.kind, i.path.name) for i in items)


def test_dry_run_deletes_nothing_and_lists_what_is_unused(tree, capsys):
    work, results = tree
    assert free_disk.main(["--work-dir", str(work), "--results-dir", str(results)]) == 0
    out = capsys.readouterr().out
    assert "dry run" in out and "legacy-db" in out and "not removed (add --include-legacy-db)" in out
    assert (work / "serve" / "v1").exists() and (work / "gene_corpus.sqlite").exists()


def test_apply_removes_old_versions_partials_and_expired_downloads_only(tree):
    work, results = tree
    assert free_disk.main(["--work-dir", str(work), "--results-dir", str(results), "--apply"]) == 0
    serve = work / "serve"
    assert not (serve / "v1").exists()                       # older than the two newest
    assert (serve / "v2").exists() and (serve / "v3").exists()
    assert (serve / "current").resolve().name == "v3"
    assert not (serve / "partial-abc").exists() and not (serve / ".presence.building.sqlite").exists()
    assert (serve / "manifest.sqlite").exists()
    assert not (results / ("a" * 20)).exists() and not (results / ("c" * 20)).exists()
    assert not (results / ("f" * 20)).exists()
    for kept in ("b", "d", "e"):                             # fresh finished, queued, running
        assert (results / (kept * 20)).exists()
    assert (results / ".quota" / "x.json").exists()           # the quota ledger is never removed
    assert (work / "gene_corpus.sqlite").exists()             # needs --include-legacy-db
    assert not (work / "per_gff" / "ann1" / "genes.sqlite.tmp123").exists()
    assert (work / "per_gff" / "ann1" / "genes.sqlite").exists()
    assert (work / "per_gff" / "gone").exists()               # needs --orphan-shards


def test_keep_versions_never_removes_current(tree):
    work, _ = tree
    (work / "serve" / "current").unlink()
    (work / "serve" / "current").symlink_to("v1")            # the oldest is live
    free_disk.main(["--work-dir", str(work), "--keep-versions", "1", "--apply"])
    assert (work / "serve" / "v1").exists() and (work / "serve" / "v3").exists()
    assert not (work / "serve" / "v2").exists()


def test_optional_categories_need_their_flags(tree):
    work, results = tree
    free_disk.main(["--work-dir", str(work), "--include-legacy-db", "--orphan-shards", "--apply"])
    assert not (work / "gene_corpus.sqlite").exists()
    assert not (work / "per_gff" / "gone").exists() and not (work / "per_gff" / "removed").exists()
    assert (work / "per_gff" / "ann1" / "genes.sqlite").exists()


def test_nothing_is_followed_through_symlinks(tree, tmp_path):
    work, results = tree
    outside = tmp_path / "outside"
    _file(outside / "keep.txt", 10)
    (work / "serve" / "partial-link").symlink_to(outside)
    free_disk.main(["--work-dir", str(work), "--apply"])
    assert (outside / "keep.txt").exists()


def test_missing_work_dir_is_an_error(tmp_path):
    assert free_disk.main(["--work-dir", str(tmp_path / "nope")]) == 2
