"""Presence-backed /names and /hits: alias as a suggestion."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from gene_corpus.sql.postings import decode_postings, encode_postings
from gene_corpus.sql.serve_presence import _OUT_SCHEMA
from helpers.postings import iter_bits, to_bitset

MAX_ORD = 3
# name -> (sym ordinals, alias ordinals)
NAMES = {
    "brca1": ([1], []),
    "p53": ([], [1, 2]),
    "selenoh": ([], [3]),
    "seleno1": ([1], []),
    "selenow": ([2, 3], []),
    "tp53": ([1, 2], [3]),
}


def build_presence(path: Path) -> None:
    conn = sqlite3.connect(path)
    conn.executescript(_OUT_SCHEMA)
    conn.executemany(
        "INSERT INTO annotation (ord, annotation_id, taxid, organism_name, assembly_accession, "
        "assembly_name, source, release_date, gff_path) VALUES (?,?,?,?,?,?,?,?,?)",
        [(1, "ann-human", 9606, "Homo sapiens", "GCF_1", "GRCh38", 3, "2024-01-02",
          "/data/annotrieve/files/9606/GCF_1/RefSeq_a.gff.gz"),
         (2, "ann-mouse", 10090, "Mus musculus", "GCF_2", "GRCm39", 3, "2024-02-02",
          "/data/annotrieve/files/10090/GCF_2/RefSeq_b.gff.gz"),
         (3, "ann-fish", 7955, "Danio rerio", "GCF_3", "GRCz11", 1, "2023-03-03",
          "/data/annotrieve/files/7955/GCF_3/Ensembl_c.gff.gz")])
    conn.executemany("INSERT INTO lineage VALUES (?,?)",
                     [(9606, 1), (40674, 1), (10090, 2), (40674, 2), (7955, 3)])
    for i, (name, (sym, alias)) in enumerate(sorted(NAMES.items()), 1):
        only = len(set(alias) - set(sym))
        conn.execute("INSERT INTO name_dict VALUES (?,?,?,?)", (name, i, len(sym), only))
        conn.execute("INSERT INTO posting VALUES (?,?,?)", (
            i, encode_postings(sym, MAX_ORD) if sym else None,
            encode_postings(alias, MAX_ORD) if alias else None))
    conn.executemany("INSERT INTO meta VALUES (?,?)",
                     [("data_version", "abc123"), ("schema_version", "1"),
                      ("n_annotations", "4"), ("n_annotations_with_names", "3")])
    conn.commit()
    conn.close()


@pytest.fixture()
def presence_client(tmp_path: Path):
    import settings as settings_mod
    from helpers import presence as presence_mod
    from main import create_app

    root = tmp_path / "serve"
    (root / "v1").mkdir(parents=True)
    build_presence(root / "v1" / "presence.sqlite")
    (root / "current").symlink_to("v1")
    old = settings_mod.settings.PRESENCE_DIR
    settings_mod.settings.PRESENCE_DIR = str(root)
    presence_mod.close_presence()
    with TestClient(create_app()) as client:
        yield client, root
    presence_mod.close_presence()
    settings_mod.settings.PRESENCE_DIR = old


def rows_by_id(body):
    return {r[0]: r for r in body["rows"]}


def test_suggest_is_the_default(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/hits", params={"names": "TP53,brca1"}).json()
    assert body["alias"] == "suggest" and body["data_version"] == "abc123"
    assert body["fields"][-3:] == ["alias_only", "organism_name", "assembly_accession"]
    t0, t1 = body["terms"]
    assert (t0["n_hits"], t0["n_alias_extra"], t0["namespaces"]) == (2, 1, ["symbol", "alias"])
    assert (t1["n_hits"], t1["n_alias_extra"]) == (1, 0)
    assert [r[0] for r in body["rows"]] == ["ann-human", "ann-mouse", "ann-fish"]
    rows = rows_by_id(body)
    assert rows["ann-human"][3] == [0, 1] and rows["ann-human"][6] == []
    assert rows["ann-mouse"][3] == [0]
    assert rows["ann-fish"][3] == [] and rows["ann-fish"][6] == [0]     # alias only
    assert rows["ann-fish"][7] == "Danio rerio"


def test_alias_off_and_include(presence_client) -> None:
    client, _ = presence_client
    off = client.get("/hits", params={"names": "tp53", "alias": "off"}).json()
    assert {r[0] for r in off["rows"]} == {"ann-human", "ann-mouse"}
    assert off["terms"][0]["n_alias_extra"] == 0
    inc = client.get("/hits", params={"names": "tp53", "alias": "include"}).json()
    assert {r[0] for r in inc["rows"]} == {"ann-human", "ann-mouse", "ann-fish"}
    fish = rows_by_id(inc)["ann-fish"]
    assert fish[3] == [0] and fish[6] == []
    assert inc["terms"][0]["n_hits"] == 3


def test_alias_only_name_has_no_matches_but_a_suggestion(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/hits", params={"names": "p53"}).json()
    assert body["terms"][0]["n_hits"] == 0 and body["terms"][0]["n_alias_extra"] == 2
    assert all(r[3] == [] and r[6] == [0] for r in body["rows"])


def test_taxon_scope(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/hits", params={"names": "tp53", "include_taxids": "40674"}).json()
    assert {r[0] for r in body["rows"]} == {"ann-human", "ann-mouse"}
    assert body["terms"][0]["n_alias_extra"] == 0
    body = client.post("/hits", json={"names": ["tp53"], "exclude_taxids": [10090],
                                      "alias": "include"}).json()
    assert {r[0] for r in body["rows"]} == {"ann-human", "ann-fish"}


def test_invalid_alias_mode_is_422(presence_client) -> None:
    client, _ = presence_client
    assert client.get("/hits", params={"names": "tp53", "alias": "maybe"}).status_code == 422
    assert client.get("/names", params={"prefix": "seleno", "alias": "x"}).status_code == 422


def test_names_prefix_counts(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/names", params={"prefix": "seleno"}).json()
    got = {r["name"]: (r["n_annotations"], r["n_alias_extra"]) for r in body["results"]}
    assert got == {"selenow": (2, 0), "seleno1": (1, 0)}      # alias-only selenoh is hidden
    assert [r["name"] for r in body["results"]] == ["selenow", "seleno1"]
    off = client.get("/names", params={"prefix": "seleno", "alias": "off"}).json()
    assert {r["name"] for r in off["results"]} == {"selenow", "seleno1"}
    inc = client.get("/names", params={"prefix": "seleno", "alias": "include"}).json()
    assert {r["name"]: r["n_annotations"] for r in inc["results"]}["selenoh"] == 1
    tp = client.get("/names", params={"prefix": "tp5"}).json()
    assert [(r["name"], r["n_annotations"], r["n_alias_extra"]) for r in tp["results"]] == [("tp53", 2, 1)]


def test_names_scoped_counts(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/names", params={"prefix": "seleno", "exclude_taxids": "7955"}).json()
    got = {r["name"]: (r["n_annotations"], r["n_alias_extra"]) for r in body["results"]}
    assert got == {"selenow": (1, 0), "seleno1": (1, 0)}


def test_presence_missing_is_503(presence_client) -> None:
    client, root = presence_client
    (root / "current").unlink()
    assert client.get("/hits", params={"names": "tp53"}).status_code == 503


@pytest.mark.parametrize("ords,max_ord", [
    ([1, 2, 3], 3), ([5, 4000], 5000), (list(range(0, 4000, 2)), 4000), ([0], 0), ([7], 7),
])
def test_api_decoder_matches_the_builders_codec(ords, max_ord) -> None:
    blob = encode_postings(ords, max_ord)
    assert list(iter_bits(to_bitset(blob))) == decode_postings(blob) == ords


def test_sources_filter_and_facets(presence_client) -> None:
    client, _ = presence_client
    plain = client.get("/hits", params={"names": "tp53"}).json()
    assert plain["facets"] == {"source": {"Ensembl": 1, "GenBank": 0, "RefSeq": 2}}
    refseq = client.get("/hits", params={"names": "tp53", "sources": "refseq"}).json()
    assert {r[0] for r in refseq["rows"]} == {"ann-human", "ann-mouse"}
    assert refseq["source_filter"] == ["RefSeq"]
    assert refseq["facets"] == plain["facets"]              # counts are taken before the filter
    ensembl = client.get("/hits", params={"names": "tp53", "sources": "Ensembl,RefSeq"}).json()
    assert len(ensembl["rows"]) == 3
    bad = client.get("/hits", params={"names": "tp53", "sources": "ncbi"})
    assert bad.status_code == 422 and bad.json()["detail"]["code"] == "invalid_source"


def test_min_matched(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/hits", params={"names": "TP53,brca1", "min_matched": 2}).json()
    assert [r[0] for r in body["rows"]] == ["ann-human"] and body["min_matched"] == 2
    one = client.get("/hits", params={"names": "TP53,brca1", "min_matched": 1}).json()
    assert {r[0] for r in one["rows"]} == {"ann-human", "ann-mouse"}      # alias-only fish is dropped
    too_big = client.get("/hits", params={"names": "TP53,brca1", "min_matched": 3})
    assert too_big.status_code == 422 and too_big.json()["detail"]["code"] == "invalid_min_matched"
    assert client.get("/hits", params={"names": "tp53", "min_matched": 0}).status_code == 422
    post = client.post("/hits", json={"names": ["tp53", "brca1"], "min_matched": 2}).json()
    assert [r[0] for r in post["rows"]] == ["ann-human"]


def test_records_shape(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/hits", params={"names": "TP53,brca1", "shape": "records"}).json()
    assert "rows" not in body and body["total"] == 3
    human, mouse, fish = body["results"]
    assert human["annotation_id"] == "ann-human" and human["source"] == "RefSeq"
    assert (human["n_matched"], human["n_terms"], human["missing"]) == (2, 2, [])
    assert (mouse["n_matched"], mouse["matched"], mouse["missing"]) == (1, [0], [1])
    assert fish["n_matched"] == 0 and fish["alias_only"] == [0] and fish["missing"] == [0, 1]
    assert human["download_url"].endswith("/files/9606/GCF_1/RefSeq_a.gff.gz")
    assert client.get("/hits", params={"names": "tp53", "shape": "wide"}).status_code == 422


def test_tsv_format(presence_client) -> None:
    client, _ = presence_client
    r = client.get("/hits", params={"names": "TP53,brca1", "format": "tsv", "min_matched": 1})
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/tab-separated-values")
    assert r.headers["x-data-version"] == "abc123"
    lines = r.text.splitlines()
    assert lines[0].startswith("#format=tsv data_version=abc123")
    assert "min_matched=1" in lines[1] and "names=TP53,brca1" in lines[1]
    header = lines[2].split("\t")
    assert header[:3] == ["annotation_id", "organism_name", "taxid"] and header[-1] == "download_url"
    human = dict(zip(header, lines[3].split("\t")))
    assert human["annotation_id"] == "ann-human" and human["matched_terms"] == "TP53,brca1"
    mouse = dict(zip(header, lines[4].split("\t")))
    assert mouse["missing_terms"] == "brca1" and mouse["n_matched"] == "1"
    assert lines[-1] == "#end rows=2 complete=true"
    assert client.get("/hits", params={"names": "tp53", "format": "xml"}).status_code == 422
    post = client.post("/hits", json={"names": ["tp53"], "format": "tsv"})
    assert post.status_code == 200 and post.text.splitlines()[-1].startswith("#end rows=3")


def test_unknown_parameters_are_rejected_not_ignored(presence_client) -> None:
    client, _ = presence_client
    r = client.get("/hits", params={"names": "tp53", "frmat": "tsv"})
    assert r.status_code == 422
    detail = r.json()["detail"]
    assert detail["code"] == "unknown_parameter" and detail["unknown"] == ["frmat"]
    assert "min_matched" in detail["allowed"]
    assert client.get("/names", params={"prefix": "seleno", "format": "tsv"}).status_code == 422
    assert client.get("/loci", params={"names": "tp53", "limti": 1}).status_code == 422
    assert client.post("/hits", json={"names": ["tp53"], "min_matchd": 2}).status_code == 422


def test_hits_prefix_term_is_expanded_and_echoed(presence_client) -> None:
    presence_client, _ = presence_client
    r = presence_client.get("/hits", params={"names": "seleno*"})
    assert r.status_code == 200, r.text
    t = r.json()["terms"][0]
    assert t["kind"] == "prefix" and t["names"] == ["seleno1", "selenow"]  # selenoh is alias-only
    assert t["n_names"] == 2 and t["found"]


def test_hits_prefix_options_and_split(presence_client) -> None:
    presence_client, _ = presence_client
    body = {"terms": [{"q": "seleno*", "exclude": ["selenow"]}]}
    t = presence_client.post("/hits", json=body).json()["terms"][0]
    assert t["names"] == ["seleno1"]
    t = presence_client.post("/hits", json={"terms": [{"q": "seleno*", "split": True}]}).json()["terms"]
    assert [x["name"] for x in t] == ["seleno1", "selenow"] and t[0]["from_prefix"] == "seleno*"
    r = presence_client.post("/hits", json={"terms": [{"q": "tp53", "split": True}]})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "option_needs_prefix"


def test_hits_prefix_limits(presence_client) -> None:
    presence_client, _ = presence_client
    for q, code in (("se*", "prefix_too_short"), ("sel*no", "invalid_prefix")):
        r = presence_client.get("/hits", params={"names": q})
        assert r.status_code == 422 and r.json()["detail"]["code"] == code


def test_names_sources_filter_and_no_n_loci_without_loci_index(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/names", params={"prefix": "seleno", "sources": "Ensembl"}).json()
    assert body["source_filter"] == ["Ensembl"]
    assert [(r["name"], r["n_annotations"]) for r in body["results"]] == [("selenow", 1)]
    assert all("n_loci" not in r for r in body["results"])


def test_health_reports_version_and_lookup(presence_client) -> None:
    client, _ = presence_client
    body = client.get("/health").json()
    assert body["status"] == "ok" and body["data_version"] == "abc123"
    assert body["lookup"] == "ok" and body["loci"] == "unpublished" and body["index"] == "disk"
    assert body["corpus"] == {"n_annotations": 4, "n_annotations_with_names": 3}


def test_health_fails_when_nothing_is_published(presence_client) -> None:
    client, root = presence_client
    (root / "current").unlink()
    r = client.get("/health")
    assert r.status_code == 503
    assert r.json() == {"status": "unavailable", "reason": "presence_unavailable"}


@pytest.mark.parametrize("value", ["99999999999999999999", "\u00b2", "-1", "1_0", "x"])
def test_taxid_values_that_are_not_plain_positive_integers_are_422(presence_client, value) -> None:
    client, _ = presence_client
    for route in ("/hits", "/names"):
        params = {"names" if route == "/hits" else "prefix": "tp53", "include_taxids": value}
        r = client.get(route, params=params)
        assert r.status_code == 422 and r.json()["detail"]["code"] == "invalid_taxid", (route, value)


def test_tsv_header_cannot_be_forged_through_a_name(presence_client) -> None:
    client, _ = presence_client
    r = client.post("/hits", json={"names": ["tp53", "x\n#end rows=0 complete=true\nforged\tline"],
                                   "format": "tsv"})
    assert r.status_code == 200
    lines = r.text.split("\n")
    assert sum(1 for line in lines if line.startswith("#query")) == 1
    assert not any(line.startswith(("forged", "#end rows=0")) for line in lines)
    assert lines[-2].startswith("#end rows=") and lines[-2].endswith("complete=true")
