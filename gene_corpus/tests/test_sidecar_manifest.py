"""Tests for gene_corpus.sql.sidecar and gene_corpus.sql.manifest (build plan M1)."""
from __future__ import annotations

import gzip
import os
import sqlite3
from pathlib import Path

import pytest

from gene_corpus.sql import per_gff
from gene_corpus.sql.manifest import (
    STATUS_FAILED,
    STATUS_INDEXED,
    STATUS_SKIPPED_COMMUNITY,
    classify_source,
    index_shard,
    open_manifest,
)
from gene_corpus.sql.sidecar import (
    SidecarError,
    is_stale,
    iter_genes,
    iter_names,
    read_header,
    shard_fingerprint,
    sidecar_path_for,
    write_sidecar,
)

BUILT_AT = "2026-01-01T00:00:00Z"


def make_shard(root: Path, annotation_id: str, *, source="RefSeq", genes=None, xrefs=None) -> Path:
    """Create ``root/<id>/genes.sqlite`` with the real per-annotation schema."""
    path = per_gff.shard_sqlite_path(root, annotation_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = per_gff.connect_per_gff(path)
    per_gff.init_schema(conn)
    per_gff.write_shard_meta(conn, {
        "annotation_id": annotation_id, "taxid": "9606",
        "assembly_accession": "GCF_000001405.40", "organism_name": "Homo sapiens",
        "source_database": source, "source_provider": "ncbi",
        "profile_id": "p", "gff_path": f"/gff/{annotation_id}.gff.gz",
    })
    genes = genes if genes is not None else [
        (1, "gene-TP53", "gene", "NC_000017.11", 7661779, 7687550, -1, "protein_coding", "TP53", None, None, 0),
        (2, "gene-OTHER", "gene", "NC_000001.11", 10, 20, 1, None, None, None, None, 0),
    ]
    xrefs = xrefs if xrefs is not None else [
        (1, "symbol", "tp53", None, "gff", 0),
        (1, "alias", "p53", None, "gff", 0),
        (1, "go", "GO:0001", None, "gff", 0),   # not a name namespace: must be ignored
        (2, "symbol", "other", None, "gff", 0),
    ]
    per_gff.insert_genes(conn, genes)
    per_gff.insert_xrefs(conn, xrefs)
    per_gff.create_indexes(conn)
    conn.commit()
    conn.close()
    return path


def raw_text(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


# --------------------------------------------------------------------------
# sidecar
# --------------------------------------------------------------------------
def test_roundtrip_names_genes_header(tmp_path):
    shard = make_shard(tmp_path, "a1")
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert info.path == sidecar_path_for(shard)
    assert (info.n_genes, info.n_name_rows, info.n_distinct_names) == (2, 3, 3)

    assert list(iter_names(info.path)) == [("other", "s", 2), ("p53", "a", 1), ("tp53", "s", 1)]
    genes = list(iter_genes(info.path))
    assert genes[0] == (1, "gene-TP53", "NC_000017.11", 7661779, 7687550, -1, "gene", "protein_coding", "TP53")
    assert genes[1][7] == "" and genes[1][8] == ""  # NULL biotype / primary_name

    header = read_header(info.path)
    assert header["annotation_id"] == "a1"
    assert header["source_database"] == "RefSeq"
    assert header["n_genes"] == "2"
    assert header["shard_fingerprint"] == shard_fingerprint(shard)


def test_non_name_namespaces_are_ignored(tmp_path):
    shard = make_shard(tmp_path, "a1")
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert all(name != "go:0001" for name, _k, _l in iter_names(info.path))
    assert "GO:0001".casefold() not in {n for n, _k, _l in iter_names(info.path)}


def test_output_is_deterministic(tmp_path):
    shard = make_shard(tmp_path, "a1")
    a = write_sidecar(shard, tmp_path / "x.gz", built_at=BUILT_AT)
    b = write_sidecar(shard, tmp_path / "y.gz", built_at=BUILT_AT)
    assert a.sidecar_fingerprint == b.sidecar_fingerprint
    assert (tmp_path / "x.gz").read_bytes() == (tmp_path / "y.gz").read_bytes()


def test_casefold_dedupe_and_two_genes_same_name(tmp_path):
    xrefs = [
        (1, "symbol", "TP53", None, "gff", 0),     # not casefolded in the shard (D4 case)
        (2, "symbol", "tp53", None, "gff", 0),     # same name, another gene
        (1, "symbol", "tp53", None, "gff", 1),     # exact duplicate after casefold? UNIQUE differs by accession only
    ]
    # UNIQUE(local_id, namespace, accession): the third row collides with the first only after casefold
    shard = make_shard(tmp_path, "a1", xrefs=xrefs)
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert list(iter_names(info.path)) == [("tp53", "s", 1), ("tp53", "s", 2)]
    assert info.n_casefold_changed == 1
    assert info.n_distinct_names == 1


def test_unsafe_names_are_skipped_and_counted(tmp_path):
    xrefs = [(1, "symbol", "good", None, "gff", 0), (1, "symbol", "bad\tname", None, "gff", 0)]
    shard = make_shard(tmp_path, "a1", xrefs=xrefs)
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert [n for n, _k, _l in iter_names(info.path)] == ["good"]
    assert info.n_skipped_names == 1


def test_gene_text_fields_cannot_break_rows(tmp_path):
    genes = [(1, "id\twith\ttabs", "gene", "chr1", 1, 2, 1, None, "na\nme", None, None, 0)]
    shard = make_shard(tmp_path, "a1", genes=genes, xrefs=[(1, "symbol", "x", None, "gff", 0)])
    info = write_sidecar(shard, built_at=BUILT_AT)
    (row,) = list(iter_genes(info.path))
    assert row[1] == "id with tabs" and row[8] == "na me"


def test_empty_shard(tmp_path):
    shard = make_shard(tmp_path, "a1", genes=[], xrefs=[])
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert (info.n_genes, info.n_name_rows) == (0, 0)
    assert list(iter_names(info.path)) == [] and list(iter_genes(info.path)) == []


def test_write_is_atomic_and_leaves_no_tmp(tmp_path):
    shard = make_shard(tmp_path, "a1")
    write_sidecar(shard, built_at=BUILT_AT)
    assert not list(shard.parent.glob("*.tmp"))


def test_truncated_file_is_detected(tmp_path):
    shard = make_shard(tmp_path, "a1")
    info = write_sidecar(shard, built_at=BUILT_AT)
    lines = raw_text(info.path).splitlines(keepends=True)
    cut = tmp_path / "cut.gz"
    with gzip.open(cut, "wt", encoding="utf-8") as fh:
        fh.writelines(lines[:-1])           # drop the #end trailer
    with pytest.raises(SidecarError):
        list(iter_names(cut))
    raw = info.path.read_bytes()
    broken = tmp_path / "broken.gz"
    broken.write_bytes(raw[: len(raw) // 2])  # cut the gzip stream itself
    with pytest.raises(SidecarError):
        list(iter_names(broken))


def test_trailer_count_mismatch_is_detected(tmp_path):
    shard = make_shard(tmp_path, "a1")
    info = write_sidecar(shard, built_at=BUILT_AT)
    text = raw_text(info.path).replace("#end names=3", "#end names=4")
    bad = tmp_path / "bad.gz"
    with gzip.open(bad, "wt", encoding="utf-8") as fh:
        fh.write(text)
    with pytest.raises(SidecarError):
        list(iter_names(bad))


def test_unsupported_version_rejected(tmp_path):
    shard = make_shard(tmp_path, "a1")
    info = write_sidecar(shard, built_at=BUILT_AT)
    text = raw_text(info.path).replace("#sidecar_version=1", "#sidecar_version=9")
    bad = tmp_path / "v9.gz"
    with gzip.open(bad, "wt", encoding="utf-8") as fh:
        fh.write(text)
    with pytest.raises(SidecarError):
        read_header(bad)


def test_shard_without_annotation_id_is_refused(tmp_path):
    shard = make_shard(tmp_path, "a1")
    conn = sqlite3.connect(shard)
    conn.execute("DELETE FROM meta WHERE key='annotation_id'")
    conn.commit()
    conn.close()
    with pytest.raises(SidecarError):
        write_sidecar(shard, built_at=BUILT_AT)


# --------------------------------------------------------------------------
# staleness (shards are rewritten in place, decision D5)
# --------------------------------------------------------------------------
def test_missing_sidecar_is_stale(tmp_path):
    shard = make_shard(tmp_path, "a1")
    assert is_stale(shard)


def test_fresh_sidecar_is_not_stale(tmp_path):
    shard = make_shard(tmp_path, "a1")
    write_sidecar(shard, built_at=BUILT_AT)
    assert not is_stale(shard)


def test_in_place_rewrite_makes_sidecar_stale_then_fresh_again(tmp_path):
    shard = make_shard(tmp_path, "a1")
    write_sidecar(shard, built_at=BUILT_AT)
    conn = sqlite3.connect(shard)
    conn.execute("INSERT INTO gene (local_id, source_gene_id, feature_type, seqid, start, end, strand) "
                 "VALUES (3, 'g3', 'gene', 'chr1', 1, 5, 1)")
    conn.commit()
    conn.close()
    assert is_stale(shard)
    write_sidecar(shard, built_at=BUILT_AT)
    assert not is_stale(shard)


def test_same_size_rewrite_is_still_detected_by_mtime(tmp_path):
    shard = make_shard(tmp_path, "a1")
    write_sidecar(shard, built_at=BUILT_AT)
    st = shard.stat()
    os.utime(shard, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert is_stale(shard)


def test_corrupt_sidecar_is_stale(tmp_path):
    shard = make_shard(tmp_path, "a1")
    write_sidecar(shard, built_at=BUILT_AT)
    sidecar_path_for(shard).write_bytes(b"not gzip")
    assert is_stale(shard)


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------
def test_classify_source():
    assert classify_source("RefSeq") is None
    assert classify_source(" ensembl ") is None
    assert classify_source("GenBank") is None
    assert classify_source("Community") == STATUS_SKIPPED_COMMUNITY
    assert classify_source("community") == STATUS_SKIPPED_COMMUNITY
    assert classify_source("CommunityRegistry") == STATUS_SKIPPED_COMMUNITY
    assert classify_source("Community Registry") == STATUS_SKIPPED_COMMUNITY
    assert classify_source("Something") == STATUS_FAILED
    assert classify_source(None) == STATUS_FAILED


def test_index_shard_registers_and_writes(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "a1")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    assert index_shard(conn, root, "a1") == STATUS_INDEXED
    row = conn.execute(
        "SELECT ord, status, n_genes, n_name_rows, shard_fingerprint, sidecar_fingerprint "
        "FROM annotation WHERE annotation_id='a1'").fetchone()
    assert row[0] == 1 and row[1] == STATUS_INDEXED and row[2] == 2 and row[3] == 3
    assert row[4].startswith("v1:") and len(row[5]) == 64
    assert sidecar_path_for(root / "a1" / "genes.sqlite").is_file()


def test_index_shard_is_idempotent_and_skips_fresh(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "a1")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    index_shard(conn, root, "a1")
    side = sidecar_path_for(root / "a1" / "genes.sqlite")
    before = side.stat().st_mtime_ns
    index_shard(conn, root, "a1")
    assert side.stat().st_mtime_ns == before          # not rewritten
    assert conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0] == 1
    index_shard(conn, root, "a1", force=True)          # force rewrites
    assert side.stat().st_mtime_ns >= before


def test_ordinals_are_stable_and_never_reused(tmp_path):
    root = tmp_path / "per_gff"
    for i in ("a1", "a2", "a3"):
        make_shard(root, i)
    conn = open_manifest(tmp_path / "manifest.sqlite")
    for i in ("a1", "a2", "a3"):
        index_shard(conn, root, i)
    ords = dict(conn.execute("SELECT annotation_id, ord FROM annotation"))
    assert ords == {"a1": 1, "a2": 2, "a3": 3}
    index_shard(conn, root, "a2", force=True)
    assert dict(conn.execute("SELECT annotation_id, ord FROM annotation")) == ords
    make_shard(root, "a4")
    index_shard(conn, root, "a4")
    assert conn.execute("SELECT ord FROM annotation WHERE annotation_id='a4'").fetchone()[0] == 4


def test_community_is_recorded_but_not_indexed(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "c1", source="Community")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    assert index_shard(conn, root, "c1") == STATUS_SKIPPED_COMMUNITY
    assert not sidecar_path_for(root / "c1" / "genes.sqlite").exists()
    assert conn.execute("SELECT COUNT(*) FROM annotation WHERE status='indexed'").fetchone()[0] == 0


def test_unknown_source_fails_visibly(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "u1", source="Mystery")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    assert index_shard(conn, root, "u1") == STATUS_FAILED
    err = conn.execute("SELECT error FROM annotation WHERE annotation_id='u1'").fetchone()[0]
    assert "Mystery" in err


def test_unreadable_shard_does_not_stop_the_run(tmp_path):
    root = tmp_path / "per_gff"
    (root / "bad").mkdir(parents=True)
    (root / "bad" / "genes.sqlite").write_bytes(b"garbage")
    make_shard(root, "ok")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    assert index_shard(conn, root, "bad") == STATUS_FAILED
    assert index_shard(conn, root, "ok") == STATUS_INDEXED
    assert index_shard(conn, root, "bad") == STATUS_FAILED   # re-run keeps one row
    assert conn.execute("SELECT COUNT(*) FROM annotation WHERE annotation_id='bad'").fetchone()[0] == 1


def test_in_place_shard_rewrite_is_picked_up_by_reindex(tmp_path):
    root = tmp_path / "per_gff"
    shard = make_shard(root, "a1")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    index_shard(conn, root, "a1")
    old = conn.execute("SELECT sidecar_fingerprint, n_genes FROM annotation").fetchone()
    c = sqlite3.connect(shard)
    c.execute("INSERT INTO gene (local_id, source_gene_id, feature_type, seqid, start, end, strand) "
              "VALUES (9, 'g9', 'gene', 'chr1', 1, 5, 1)")
    c.execute("INSERT INTO gene_xref VALUES (9, 'symbol', 'newname', NULL, 'gff', 0)")
    c.commit()
    c.close()
    index_shard(conn, root, "a1")
    new = conn.execute("SELECT sidecar_fingerprint, n_genes FROM annotation").fetchone()
    assert new[0] != old[0] and new[1] == old[1] + 1
    assert ("newname", "s", 9) in set(iter_names(sidecar_path_for(shard)))


def test_large_shard_uses_file_scratch_and_cleans_up(tmp_path, monkeypatch):
    from gene_corpus.sql import sidecar as sc

    monkeypatch.setattr(sc, "_MEMORY_ROWS", 1)       # force the file-backed scratch table
    monkeypatch.setattr(sc, "_BATCH", 2)
    shard = make_shard(tmp_path / "per_gff", "a1")
    tmpd = tmp_path / "scratch"
    info = write_sidecar(shard, built_at=BUILT_AT, tmp_dir=tmpd)
    assert list(iter_names(info.path)) == [("other", "s", 2), ("p53", "a", 1), ("tp53", "s", 1)]
    assert not list(tmpd.glob("*")), "scratch file must be removed"
    ref = write_sidecar(shard, tmp_path / "ref.gz", built_at=BUILT_AT)
    assert ref.sidecar_fingerprint == info.sidecar_fingerprint   # same bytes either way


def test_failed_write_leaves_no_tmp_and_keeps_old_sidecar(tmp_path):
    shard = make_shard(tmp_path, "a1")
    good = write_sidecar(shard, built_at=BUILT_AT)
    conn = sqlite3.connect(shard)
    conn.execute("DELETE FROM meta WHERE key='annotation_id'")
    conn.commit()
    conn.close()
    with pytest.raises(SidecarError):
        write_sidecar(shard, built_at=BUILT_AT)
    assert not list(shard.parent.glob("*.tmp"))
    assert good.path.is_file()


def test_genes_section_holds_only_named_genes(tmp_path):
    genes = [(g, f"g{g}", "gene", "chr1", g, g + 5, 1, None, f"N{g}", None, None, 0)
             for g in range(1, 1201)]
    xrefs = [(g, "symbol", f"name{g}", None, "gff", 0) for g in range(1, 1201, 3)]   # every 3rd gene
    xrefs.append((7, "alias", "seven", None, "gff", 0))                              # gene 7 has two names
    shard = make_shard(tmp_path, "a1", genes=genes, xrefs=xrefs)
    info = write_sidecar(shard, built_at=BUILT_AT)
    named = sorted({lid for _n, _k, lid in iter_names(info.path)})
    assert [g[0] for g in iter_genes(info.path)] == named
    assert info.n_genes == 1200 and info.n_genes_named == len(named) == 400
    assert read_header(info.path)["n_genes_named"] == "400"


def test_shard_without_any_names_writes_no_genes(tmp_path):
    genes = [(g, f"g{g}", "gene", "chr1", g, g + 5, 1, None, None, None, None, 0) for g in range(1, 50)]
    shard = make_shard(tmp_path, "a1", genes=genes, xrefs=[(1, "go", "GO:1", None, "gff", 0)])
    info = write_sidecar(shard, built_at=BUILT_AT)
    assert (info.n_genes, info.n_genes_named, info.n_name_rows) == (49, 0, 0)
    assert list(iter_genes(info.path)) == []


def test_community_registry_shard_is_never_indexed(tmp_path):
    root = tmp_path / "per_gff"
    make_shard(root, "c1", source="CommunityRegistry")
    conn = open_manifest(tmp_path / "manifest.sqlite")
    assert index_shard(conn, root, "c1") == STATUS_SKIPPED_COMMUNITY
    assert not sidecar_path_for(root / "c1" / "genes.sqlite").exists()
