"""Presence-backed /names and /hits: alias as a suggestion (API-DESIGN section 3)."""
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
        [(1, "ann-human", 9606, "Homo sapiens", "GCF_1", "GRCh38", 3, "2024-01-02", None),
         (2, "ann-mouse", 10090, "Mus musculus", "GCF_2", "GRCm39", 3, "2024-02-02", None),
         (3, "ann-fish", 7955, "Danio rerio", "GCF_3", "GRCz11", 1, "2023-03-03", None)])
    conn.executemany("INSERT INTO lineage VALUES (?,?)",
                     [(9606, 1), (40674, 1), (10090, 2), (40674, 2), (7955, 3)])
    for i, (name, (sym, alias)) in enumerate(sorted(NAMES.items()), 1):
        only = len(set(alias) - set(sym))
        conn.execute("INSERT INTO name_dict VALUES (?,?,?,?)", (name, i, len(sym), only))
        conn.execute("INSERT INTO posting VALUES (?,?,?)", (
            i, encode_postings(sym, MAX_ORD) if sym else None,
            encode_postings(alias, MAX_ORD) if alias else None))
    conn.executemany("INSERT INTO meta VALUES (?,?)",
                     [("data_version", "abc123"), ("schema_version", "1")])
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
    assert got == {"selenow": (2, 0), "seleno1": (1, 0), "selenoh": (0, 1)}
    assert [r["name"] for r in body["results"]] == ["selenow", "seleno1", "selenoh"]
    off = client.get("/names", params={"prefix": "seleno", "alias": "off"}).json()
    assert {r["name"] for r in off["results"]} == {"selenow", "seleno1"}
    inc = client.get("/names", params={"prefix": "seleno", "alias": "include"}).json()
    assert {r["name"]: r["n_annotations"] for r in inc["results"]}["selenoh"] == 1


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
