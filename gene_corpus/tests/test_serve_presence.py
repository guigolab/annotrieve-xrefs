"""Tests for gene_corpus.sql.postings / serve_presence / serve_build CLI."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.sql.postings import (
    TAG_BITSET,
    TAG_VARINT,
    bitset_nbytes,
    decode_bitset,
    decode_postings,
    encode_postings,
)
from gene_corpus.sql.serve_presence import BuildError, build_presence, publish_version


# --------------------------------------------------------------------------
# codec
# --------------------------------------------------------------------------
def test_roundtrip_small_uses_varint():
    ords = [1, 5, 300, 20000]
    blob = encode_postings(ords, 20000)
    assert blob[0] == TAG_VARINT
    assert decode_postings(blob) == ords
    assert decode_bitset(blob) == sum(1 << o for o in ords)


def test_roundtrip_dense_uses_bitset():
    ords = list(range(1, 3000, 2))
    blob = encode_postings(ords, 3000)
    assert blob[0] == TAG_BITSET
    assert len(blob) == 1 + bitset_nbytes(3000)
    assert decode_postings(blob) == ords
    assert decode_bitset(blob) == sum(1 << o for o in ords)


def test_chooses_smaller_encoding():
    # 5 entries spread widely: varint is 3 bytes each = 15 > nbytes(1000)=126? no -> varint
    blob = encode_postings([10, 400, 700, 900, 1000], 1000)
    assert blob[0] == TAG_VARINT
    assert len(blob) <= 1 + bitset_nbytes(1000)


def test_single_and_first_ordinal_zero():
    assert decode_postings(encode_postings([0], 10)) == [0]
    assert decode_postings(encode_postings([7], 7)) == [7]


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        encode_postings([], 10)
    with pytest.raises(ValueError):
        encode_postings([3, 3], 10)
    with pytest.raises(ValueError):
        encode_postings([5, 2], 10)
    with pytest.raises(ValueError):
        encode_postings([11], 10)
    with pytest.raises(ValueError):
        decode_postings(bytes([0x7F, 1]))


# --------------------------------------------------------------------------
# fixture source DB
# --------------------------------------------------------------------------
def _make_source(path: Path, *, xrefs: list[tuple[str, str, int]]) -> None:
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY, annotation_id TEXT UNIQUE NOT NULL,
            taxid INTEGER NOT NULL, assembly_accession TEXT NOT NULL, organism_name TEXT,
            source_database TEXT NOT NULL, source_provider TEXT, profile_id TEXT NOT NULL,
            gff_path TEXT, release_date TEXT, assembly_name TEXT);
        CREATE TABLE annotation_lineage (
            annotation_key INTEGER NOT NULL, taxid INTEGER NOT NULL,
            PRIMARY KEY (annotation_key, taxid)) WITHOUT ROWID;
        CREATE TABLE xref_file (
            namespace TEXT NOT NULL, accession TEXT NOT NULL, annotation_key INTEGER NOT NULL,
            PRIMARY KEY (namespace, accession, annotation_key)) WITHOUT ROWID;
        """
    )
    rows = [
        (1, "aaa", 9606, "GCF_1", "Homo sapiens", "RefSeq", None, "p", "a.gff.gz", "2024-01-01", "GRCh38"),
        (2, "bbb", 10090, "GCF_2", "Mus musculus", "Ensembl", None, "p", "b.gff.gz", "2024-02-01", "GRCm39"),
        (3, "ccc", 7955, "GCA_3", "Danio rerio", "GenBank", None, "p", "c.gff.gz", None, None),
    ]
    conn.executemany("INSERT INTO annotation VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
    conn.executemany(
        "INSERT INTO annotation_lineage VALUES (?,?)",
        [(1, 9606), (1, 40674), (2, 10090), (2, 40674), (3, 7955), (99, 1)],  # 99 = orphan
    )
    conn.executemany("INSERT INTO xref_file VALUES (?,?,?)", xrefs)
    conn.commit()
    conn.close()


XREFS = [
    ("symbol", "brca1", 1), ("symbol", "brca1", 2), ("symbol", "brca1", 3),
    ("symbol", "tp53", 1), ("symbol", "tp53", 2),
    ("alias", "tp53", 3),            # tp53 symbol in 1,2; alias-only in 3
    ("alias", "p53", 1), ("alias", "p53", 2),   # alias-only name
    ("symbol", "orphanonly", 42),    # annotation_key 42 is not in annotation
    ("symbol", "mixed", 2), ("symbol", "mixed", 42),
    ("symbol", "tp53b", 3),
]


@pytest.fixture()
def source(tmp_path: Path) -> Path:
    p = tmp_path / "gene_corpus.sqlite"
    _make_source(p, xrefs=XREFS)
    return p


def _out(tmp_path: Path) -> Path:
    return tmp_path / "presence.sqlite"


def _names(out: Path) -> dict[str, tuple]:
    conn = sqlite3.connect(out)
    try:
        return {
            r[0]: tuple(r[1:])
            for r in conn.execute(
                "SELECT name, name_id, n_sym, n_alias_only FROM name_dict ORDER BY name")
        }
    finally:
        conn.close()


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------
def test_build_contents(source, tmp_path):
    out = _out(tmp_path)
    stats = build_presence(source, out, log=lambda m: None)
    names = _names(out)
    # orphan-only name is dropped; mixed keeps only its valid annotation
    assert set(names) == {"brca1", "mixed", "p53", "tp53", "tp53b"}
    assert names["tp53"][1:] == (2, 1)        # n_sym=2, n_alias_only=1 (annotation 3)
    assert names["p53"][1:] == (0, 2)         # alias-only name
    assert names["mixed"][1:] == (1, 0)
    assert names["brca1"][1:] == (3, 0)
    # name_ids follow sorted name order
    assert [v[0] for _k, v in sorted(names.items())] == [1, 2, 3, 4, 5]
    assert stats["orphan_rows_dropped"] == 2
    assert stats["n_annotations"] == 3 and stats["partial"] is False

    conn = sqlite3.connect(out)
    try:
        sym, alias = conn.execute(
            "SELECT sym, alias FROM posting WHERE name_id = ?", (names["tp53"][0],)).fetchone()
        assert decode_postings(sym) == [1, 2]
        assert decode_postings(alias) == [3]
        assert conn.execute(
            "SELECT sym FROM posting WHERE name_id = ?", (names["p53"][0],)).fetchone()[0] is None
        meta = dict(conn.execute("SELECT key, value FROM meta"))
        assert meta["partial"] == "0" and meta["n_annotations"] == "3"
        assert meta["max_ord"] == "3" and len(meta["data_version"]) == 16
        ann = conn.execute(
            "SELECT ord, annotation_id, source, assembly_name FROM annotation ORDER BY ord"
        ).fetchall()
        assert ann == [(1, "aaa", 3, "GRCh38"), (2, "bbb", 1, "GRCm39"), (3, "ccc", 2, None)]
        lin = conn.execute("SELECT taxid, ord FROM lineage ORDER BY taxid, ord").fetchall()
        assert (1, 99) not in lin and (40674, 1) in lin and len(lin) == 5
    finally:
        conn.close()


def test_prefix_build_is_partial(source, tmp_path):
    out = _out(tmp_path)
    stats = build_presence(source, out, name_prefix="TP5", log=lambda m: None)
    assert stats["partial"] is True
    assert set(_names(out)) == {"tp53", "tp53b"}
    conn = sqlite3.connect(out)
    try:
        assert conn.execute("SELECT value FROM meta WHERE key='partial'").fetchone()[0] == "1"
    finally:
        conn.close()


def test_limit_names(source, tmp_path):
    out = _out(tmp_path)
    stats = build_presence(source, out, limit_names=2, log=lambda m: None)
    assert stats["n_names"] == 2 and stats["partial"] is True
    assert list(_names(out)) == ["brca1", "mixed"]


def test_dense_name_uses_bitset(tmp_path):
    src = tmp_path / "gene_corpus.sqlite"
    conn = sqlite3.connect(src)
    conn.executescript(
        """
        CREATE TABLE annotation (annotation_key INTEGER PRIMARY KEY, annotation_id TEXT,
            taxid INT, assembly_accession TEXT, organism_name TEXT, source_database TEXT,
            source_provider TEXT, profile_id TEXT, gff_path TEXT, release_date TEXT, assembly_name TEXT);
        CREATE TABLE annotation_lineage (annotation_key INT, taxid INT);
        CREATE TABLE xref_file (namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)) WITHOUT ROWID;
        """
    )
    n = 400
    conn.executemany("INSERT INTO annotation VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                     [(k, f"id{k}", 9606, "A", "x", "RefSeq", None, "p", None, None, None)
                      for k in range(1, n + 1)])
    conn.executemany("INSERT INTO xref_file VALUES ('symbol','actb',?)",
                     [(k,) for k in range(1, n + 1)])
    conn.executemany("INSERT INTO xref_file VALUES ('symbol','rare',?)", [(5,), (399,)])
    conn.commit()
    conn.close()
    out = _out(tmp_path)
    stats = build_presence(src, out, log=lambda m: None)
    assert stats["bitset_blobs"] == 1 and stats["varint_blobs"] == 1
    c = sqlite3.connect(out)
    try:
        blob = c.execute(
            "SELECT sym FROM posting p JOIN name_dict d USING(name_id) WHERE d.name='actb'"
        ).fetchone()[0]
        assert blob[0] == TAG_BITSET
        assert decode_postings(blob) == list(range(1, n + 1))
    finally:
        c.close()


def test_unsorted_stream_is_refused(source, tmp_path, monkeypatch):
    import gene_corpus.sql.serve_presence as sp

    real = sp.groupby

    def bad_groupby(it, key=None):
        groups = list(real(it, key=key))
        groups.reverse()
        return iter(groups)

    monkeypatch.setattr(sp, "groupby", bad_groupby)
    with pytest.raises(BuildError):
        build_presence(source, _out(tmp_path), log=lambda m: None)


def test_verification_detects_tampering(source, tmp_path, monkeypatch):
    import gene_corpus.sql.serve_presence as sp

    real_encode = sp.encode_postings

    def wrong_encode(ords, max_ord):
        if len(ords) == 2:  # corrupt tp53 symbol postings
            return real_encode([ords[0]], max_ord)
        return real_encode(ords, max_ord)

    monkeypatch.setattr(sp, "encode_postings", wrong_encode)
    with pytest.raises(BuildError):
        build_presence(source, _out(tmp_path), verify_sample=100, log=lambda m: None)


def test_no_annotations_is_error(tmp_path):
    src = tmp_path / "gene_corpus.sqlite"
    _make_source(src, xrefs=[])
    c = sqlite3.connect(src)
    c.execute("DELETE FROM annotation")
    c.commit()
    c.close()
    with pytest.raises(BuildError):
        build_presence(src, _out(tmp_path), log=lambda m: None)


# --------------------------------------------------------------------------
# publish / CLI
# --------------------------------------------------------------------------
def test_publish_switches_current_symlink(source, tmp_path):
    out_dir = tmp_path / "serve"
    out_dir.mkdir()
    built = out_dir / ".b.sqlite"
    stats = build_presence(source, built, log=lambda m: None)
    final = publish_version(out_dir, built, stats)
    assert final.is_file() and not built.exists()
    current = out_dir / "current"
    assert current.is_symlink()
    assert (current / "presence.sqlite").is_file()
    assert json.loads((current / "build_report.json").read_text())["n_names"] == 5
    # a second publish replaces the link atomically and keeps the old directory
    built2 = out_dir / ".b2.sqlite"
    stats2 = dict(build_presence(source, built2, name_prefix="b", log=lambda m: None))
    stats2["data_version"] = "ffffffffffffffff"
    publish_version(out_dir, built2, stats2)
    assert current.resolve().name == "vffffffffffffffff"
    assert final.is_file()


def test_cli_full_and_partial(source, tmp_path):
    work = source.parent
    assert serve_build_main(["--work-dir", str(work)]) == 0
    serve = work / "serve"
    assert (serve / "current" / "presence.sqlite").is_file()
    assert serve_build_main(["--work-dir", str(work), "--name-prefix", "tp",
                             "--out-dir", str(tmp_path / "trial")]) == 0
    trial = tmp_path / "trial"
    assert not (trial / "current").exists()
    assert list(trial.glob("partial-*/presence.sqlite"))


def test_cli_missing_source(tmp_path):
    assert serve_build_main(["--work-dir", str(tmp_path)]) == 2


def test_verify_checks_every_name_when_sample_exceeds_table(source, tmp_path):
    stats = build_presence(source, _out(tmp_path), verify_sample=1000, log=lambda m: None)
    assert stats["verified_names"] == stats["n_names"] == 5
