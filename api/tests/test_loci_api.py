"""/loci over a real sidecar-built presence + loci index."""
from __future__ import annotations

from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest
from fastapi.testclient import TestClient

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.tests.test_presence_from_sidecars import corpus


@pytest.fixture()
def loci_client(tmp_path: Path):
    import settings as settings_mod
    from helpers import presence as presence_mod
    from main import create_app

    work = corpus(tmp_path)
    assert serve_build_main(["build", "--work-dir", str(work)]) == 0
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 0
    old = settings_mod.settings.PRESENCE_DIR
    settings_mod.settings.PRESENCE_DIR = str(work / "serve")
    presence_mod.close_presence()
    with TestClient(create_app()) as client:
        yield client, work
    presence_mod.close_presence()
    settings_mod.settings.PRESENCE_DIR = old


def test_json_rows_and_dry_run(loci_client) -> None:
    client, _ = loci_client
    body = client.get("/loci", params={"names": "TP53", "alias": "off"}).json()
    assert len(body["rows"]) == 5 and body["next_cursor"] is None
    r = body["rows"][0]
    assert (r["seqid"], r["start"], r["end"], r["strand"]) == ("chr1", 100, 150, "+")
    assert r["symbol"] == "tp53" and r["matched_via"] == "symbol" and r["terms"] == [0]
    assert r["locus_ref"] == f"{r['annotation_id']}:1" and r["gene_id"] == "gene-1"
    q = parse_qs(urlparse(r["gff_url"]).query)
    assert q == {"region": ["chr1"], "start": ["99"], "end": ["150"]}      # 0-based half-open
    dry = client.get("/loci", params={"names": "tp53", "dry_run": "true"}).json()
    assert dry["n_annotations"] == 5 and dry["n_loci_upper_bound"] == 7 and "rows" not in dry


def test_alias_modes_and_flank(loci_client) -> None:
    client, _ = loci_client
    sug = client.get("/loci", params={"names": "tp53", "flank_bp": 10}).json()
    assert len(sug["rows"]) == 7                       # 5 symbol genes + 2 alias-only genes
    via = {(r["local_id"], r["matched_via"]) for r in sug["rows"]}
    assert via == {(1, "symbol"), (2, "alias")}
    sym = next(r for r in sug["rows"] if r["local_id"] == 1)
    assert parse_qs(urlparse(sym["gff_url"]).query)["start"] == ["89"]
    assert parse_qs(urlparse(sym["gff_url"]).query)["end"] == ["160"]
    assert len(client.get("/loci", params={"names": "tp53", "alias": "off"}).json()["rows"]) == 5


def test_one_row_per_gene_with_several_terms_and_order(loci_client) -> None:
    client, _ = loci_client
    body = client.get("/loci", params={"names": "tp53,p53"}).json()
    gene2 = [r for r in body["rows"] if r["local_id"] == 2]
    assert gene2 and all(r["terms"] == [0, 1] or r["terms"] == [0] or r["terms"] == [1] for r in gene2)
    assert len({(r["annotation_id"], r["local_id"]) for r in body["rows"]}) == len(body["rows"])
    ids = [r["annotation_id"] for r in body["rows"]]
    assert ids == sorted(ids)                           # ordinals follow sorted shard order


def test_pagination_cursor_walks_every_row_once(loci_client) -> None:
    client, _ = loci_client
    seen, cursor = [], None
    for _ in range(20):
        params = {"names": "tp53", "limit": 3}
        if cursor:
            params["cursor"] = cursor
        body = client.get("/loci", params=params).json()
        seen += [r["locus_ref"] for r in body["rows"]]
        cursor = body["next_cursor"]
        if not cursor:
            break
    full = [r["locus_ref"] for r in client.get("/loci", params={"names": "tp53", "limit": 200}).json()["rows"]]
    assert seen == full and len(set(seen)) == 7
    assert client.get("/loci", params={"names": "tp53", "cursor": "bad"}).status_code == 422
    c1 = client.get("/loci", params={"names": "tp53", "limit": 1}).json()["next_cursor"]
    assert client.get("/loci", params={"names": "brca1", "cursor": c1}).status_code == 422


def test_tsv_ndjson_streams_with_trailer(loci_client) -> None:
    client, _ = loci_client
    tsv = client.get("/loci", params={"names": "tp53", "format": "tsv", "alias": "off"})
    lines = tsv.text.splitlines()
    assert lines[0].startswith("#format=tsv data_version=") and lines[1].startswith("#query names=tp53")
    assert lines[2].split("\t")[0] == "annotation_id" and lines[-1] == "#end rows=5 complete=true"
    assert len([ln for ln in lines if not ln.startswith("#") and not ln.startswith("annotation_id")]) == 5
    bad = client.get("/loci", params={"names": "tp53", "format": "bed"})
    assert bad.status_code == 422
    nd = client.get("/loci", params={"names": "tp53", "format": "ndjson"}).text.splitlines()
    assert len(nd) == 8 and '"end"' in nd[-1] and '"complete": true' in nd[-1]
    cut = client.get("/loci", params={"names": "tp53", "format": "tsv", "limit": 2}).text.splitlines()
    assert cut[-1] == "#end rows=2 complete=false"


def test_scope_filters_and_errors(loci_client) -> None:
    client, _ = loci_client
    ids = [r["annotation_id"] for r in client.get("/loci", params={"names": "tp53", "alias": "off"}).json()["rows"]]
    one = client.get("/loci", params={"names": "tp53", "annotation_ids": ids[1]}).json()
    assert {r["annotation_id"] for r in one["rows"]} == {ids[1]}
    assert client.get("/loci", params={"names": "tp53", "include_taxids": "1"}).json()["rows"] != []
    assert client.get("/loci", params={"names": "tp53", "exclude_taxids": "9606"}).json()["rows"] == []
    assert client.get("/loci", params={"names": "nosuchgene"}).json()["terms"][0]["found"] is False
    assert client.get("/loci", params={"names": "tp53", "format": "xml"}).status_code == 422
    assert client.get("/loci", params={"names": "tp53", "flank_bp": -1}).status_code == 422
    assert client.get("/loci", params={"names": "tp53", "alias": "x"}).status_code == 422
    post = client.post("/loci", json={"names": ["tp53"], "alias": "off", "format": "ndjson"})
    assert post.status_code == 200 and post.text.count("\n") == 6


def test_case_sensitive_uses_the_primary_name_spelling(loci_client) -> None:
    client, work = loci_client
    # shard primary names are None in the fixture genes, so no row has TP53 spelled exactly
    body = client.get("/loci", params={"names": "tp53", "case_sensitive": "true"}).json()
    assert body["rows"] == []


def test_export_cap_returns_413_for_unscoped_queries(loci_client, monkeypatch) -> None:
    client, _ = loci_client
    import settings as settings_mod
    monkeypatch.setattr(settings_mod.settings, "LOCI_MAX_ROWS", 3)
    r = client.get("/loci", params={"names": "tp53", "format": "tsv"})
    assert r.status_code == 413 and r.json()["detail"]["code"] == "export_too_large"
    assert client.get("/loci", params={"names": "tp53", "format": "tsv", "limit": 10}).status_code == 200


def test_missing_loci_file_is_503(loci_client) -> None:
    client, work = loci_client
    (work / "serve" / "current" / "loci.sqlite").unlink()
    assert client.get("/loci", params={"names": "tp53"}).status_code == 503


BAD_LOCI_PARAMS = [
    {"names": "tp53", "include_taxids": "abc"},
    {"names": "tp53", "exclude_taxids": "1,,x"},
    {"names": "tp53", "alias": "bogus"},
    {"names": "tp53", "flank_bp": -5},
    {"names": "tp53", "flank_bp": 10**12},
    {"names": "tp53", "format": "xml"},
    {"names": "tp53", "limit": -1},
    {"names": "tp53", "limit": 0},
    {"names": "tp53", "cursor": "not-a-cursor"},
    {"names": "tp53", "cursor": "%%%"},
    {"names": "tp5*"},
    {"names": ""},
    {"names": ",,,"},
    {"names": "tp53", "annotation_ids": ""},
    {"names": "x" * 5000},
]


@pytest.mark.parametrize("params", BAD_LOCI_PARAMS)
def test_bad_loci_parameters_never_cause_a_server_error(loci_client, params) -> None:
    client, _ = loci_client
    r = client.get("/loci", params=params)
    assert r.status_code < 500, (params, r.status_code, r.text[:300])


@pytest.mark.parametrize("params", [
    {"prefix": "tp"}, {"prefix": "tp5*"}, {"prefix": "tp5", "include_taxids": "x"},
    {"prefix": "tp5", "alias": "zzz"}, {"prefix": "%"}, {"prefix": "_" * 3}])
def test_bad_names_parameters_never_cause_a_server_error(loci_client, params) -> None:
    client, _ = loci_client
    assert client.get("/names", params=params).status_code < 500


@pytest.mark.parametrize("params", [
    {"names": "tp53", "include_taxids": "x"}, {"names": "tp53", "alias": "zzz"},
    {"names": ",,"}, {"names": "tp5*"}, {"names": "a," * 5000}])
def test_bad_hits_parameters_never_cause_a_server_error(loci_client, params) -> None:
    client, _ = loci_client
    assert client.get("/hits", params=params).status_code < 500
    assert client.post("/hits", json={"names": params["names"].split(",")}).status_code < 500


# ---- the legacy /annotations routes must agree with the new index ---------------------------

@pytest.fixture()
def shard_client(loci_client):
    import settings as settings_mod
    client, work = loci_client
    old = settings_mod.settings.GENE_CORPUS_DIR
    settings_mod.settings.GENE_CORPUS_DIR = str(work)
    try:
        yield client, work
    finally:
        settings_mod.settings.GENE_CORPUS_DIR = old


def test_every_locus_ref_resolves_to_the_same_gene_in_the_annotation_routes(shard_client) -> None:
    client, _ = shard_client
    rows = client.get("/loci", params={"names": "tp53", "alias": "include"}).json()["rows"]
    assert rows
    for r in rows:
        ann, local = r["locus_ref"].rsplit(":", 1)
        gene = client.get(f"/annotations/{ann}/genes/{local}")
        assert gene.status_code == 200, (r["locus_ref"], gene.text)
        g = gene.json()
        assert (g["local_id"], g["source_gene_id"], g["seqid"], g["start"], g["end"]) == (
            int(local), r["gene_id"], r["seqid"], r["start"], r["end"]), r["locus_ref"]
        assert g["strand"] == r["strand"] and g["locus_ref"] == r["locus_ref"]
        assert client.get(f"/annotations/{ann}/genes/{local}/xrefs").status_code == 200


def test_annotation_genes_filter_finds_the_genes_loci_lists(shard_client) -> None:
    client, _ = shard_client
    rows = client.get("/loci", params={"names": "tp53", "alias": "include"}).json()["rows"]
    by_ann: dict[str, set[int]] = {}
    for r in rows:
        ann, local = r["locus_ref"].rsplit(":", 1)
        by_ann.setdefault(ann, set()).add(int(local))
    for ann, expected in by_ann.items():
        body = client.get(f"/annotations/{ann}/genes", params={"names": "tp53"}).json()
        assert {g["local_id"] for g in body["results"]} == expected, ann
        assert all(g["locus_ref"] == f"{ann}:{g['local_id']}" for g in body["results"])


@pytest.mark.parametrize("path,params", [
    ("/annotations/%s/genes", {"limit": 0}), ("/annotations/%s/genes", {"limit": 10**9}),
    ("/annotations/%s/genes", {"offset": -1}), ("/annotations/%s/genes", {"names": "a*"}),
    ("/annotations/%s/genes", {"names": "a", "curies": "x:1"}),
    ("/annotations/%s/genes", {"curies": "nocolon"}), ("/annotations/%s/genes", {"match": "zzz"}),
    ("/annotations/%s/genes/-1", {}), ("/annotations/%s/genes/abc", {}),
    ("/annotations/%s/genes/99999999999999999999", {}),
    ("/annotations/%s/genes/1/xrefs", {"offset": 10**12}),
    ("/annotations/%s/namespaces/nope/accessions", {}),
    ("/annotations/%s/namespaces/a%%2Fb/accessions", {"q": "x"}),
])
def test_annotation_routes_never_cause_a_server_error(shard_client, path, params) -> None:
    client, _ = shard_client
    ann = client.get("/loci", params={"names": "tp53"}).json()["rows"][0]["annotation_id"]
    assert client.get(path % ann, params=params).status_code < 500


@pytest.mark.parametrize("ann", ["nope", "..", "a" * 500, "%00", "..%2F..%2Fetc"])
def test_unknown_or_hostile_annotation_ids_are_404_not_500(shard_client, ann) -> None:
    client, _ = shard_client
    for path in ("namespaces", "genes", "genes/1", "genes/1/xrefs"):
        assert client.get(f"/annotations/{ann}/{path}").status_code in (404, 422)


def test_loci_prefix_term_expands_and_reports_matched_names(loci_client) -> None:
    client, _ = loci_client
    body = client.get("/loci", params={"names": "uniq*", "alias": "off"}).json()
    t = body["terms"][0]
    assert t["kind"] == "prefix" and t["names"] == [f"uniq{i}" for i in range(5)]
    assert len(body["rows"]) == 5
    assert {r["symbol"] for r in body["rows"]} == set(t["names"])
    assert all(r["matched_names"] == [r["symbol"]] and r["terms"] == [0] for r in body["rows"])
    dry = client.get("/loci", params={"names": "uniq*", "dry_run": "true"}).json()
    assert dry["n_annotations"] == 5 and "rows" not in dry


def test_loci_split_prefix_and_pagination(loci_client) -> None:
    client, _ = loci_client
    body = client.post("/loci", json={"terms": [{"q": "uniq*", "split": True, "exclude": ["uniq4"]}],
                                      "alias": "off", "limit": 2}).json()
    assert [t["name"] for t in body["terms"]] == [f"uniq{i}" for i in range(4)]
    seen, cur = [], body
    while True:
        seen += [r["locus_ref"] for r in cur["rows"]]
        if not cur["next_cursor"]:
            break
        cur = client.post("/loci", json={"terms": [{"q": "uniq*", "split": True,
                          "exclude": ["uniq4"]}], "alias": "off", "limit": 2,
                          "cursor": cur["next_cursor"]}).json()
    assert len(seen) == 4 == len(set(seen))


def test_loci_prefix_errors_and_case_sensitive(loci_client) -> None:
    client, _ = loci_client
    r = client.get("/loci", params={"names": "un*"})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "prefix_too_short"
    rows = client.get("/loci", params={"names": "UNIQ*", "case_sensitive": "true"}).json()["rows"]
    assert rows == []
    # fixture genes have no primary_name, so a case-sensitive prefix cannot match their spelling
    assert client.get("/loci", params={"names": "uniq*", "case_sensitive": "true",
                                       "alias": "off"}).json()["rows"] == []


def test_names_n_loci_and_sources_filter(loci_client) -> None:
    client, _ = loci_client
    body = client.get("/names", params={"prefix": "tp5"}).json()
    row = body["results"][0]
    assert row["name"] == "tp53" and row["n_loci"] == 5            # symbol gene rows only
    inc = client.get("/names", params={"prefix": "tp5", "alias": "include"}).json()["results"][0]
    assert inc["n_loci"] == 7
    r = client.get("/names", params={"prefix": "tp5", "sources": "Nope"})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "invalid_source"
    scoped = client.get("/names", params={"prefix": "tp5", "sources": "RefSeq,Ensembl,GenBank"}
                        ).json()["results"][0]
    assert scoped["n_loci"] == 5 and scoped["n_annotations"] == row["n_annotations"]


def test_a_bad_cursor_is_a_422_with_a_code(loci_client) -> None:
    client = loci_client[0]
    r = client.get("/loci", params={"names": "tp53", "cursor": "not-a-cursor"})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "invalid_cursor"


def test_loci_sources_filter_keeps_only_those_annotations(loci_client) -> None:
    client, _ = loci_client
    allrows = client.get("/loci", params={"names": "tp53", "limit": 200}).json()
    ids_all = {r["annotation_id"] for r in allrows["rows"]}
    facets = client.get("/hits", params={"names": "tp53"}).json()["facets"]["source"]
    seen = set()
    for label, n in facets.items():
        want = {r[0] for r in client.get("/hits", params={"names": "tp53", "sources": label}).json()["rows"]}
        got = client.get("/loci", params={"names": "tp53", "sources": label, "limit": 200}).json()
        assert got["source_filter"] == [label]
        assert {r["annotation_id"] for r in got["rows"]} <= want
        assert {r["annotation_id"] for r in got["rows"]} == (want & ids_all)
        seen |= {r["annotation_id"] for r in got["rows"]}
    assert seen == ids_all                                   # every source together = no filter
    both = client.get("/loci", params={"names": "tp53", "sources": ",".join(facets), "limit": 200}).json()
    assert {r["annotation_id"] for r in both["rows"]} == ids_all
    assert client.post("/loci", json={"names": ["tp53"], "sources": ["Nowhere"]}).status_code == 422
    bad = client.get("/loci", params={"names": "tp53", "sources": "Nowhere"})
    assert bad.status_code == 422 and bad.json()["detail"]["code"] == "invalid_source"


def test_a_cursor_from_one_source_filter_is_refused_under_another(loci_client) -> None:
    client, _ = loci_client
    first = client.get("/loci", params={"names": "tp53", "limit": 1}).json()
    assert first["next_cursor"]
    r = client.get("/loci", params={"names": "tp53", "limit": 1, "cursor": first["next_cursor"],
                                    "sources": "RefSeq"})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "cursor_filter_mismatch"


def test_loci_tsv_provenance_is_one_line_whatever_the_names(tmp_path) -> None:
    from services.loci_formats import _provenance
    from services.presence_loci import LociQuery

    q = LociQuery("suggest", False, 0, [], [], [], [], [{"input": "a\nb\tc"}], [], set(), "v1")
    lines = _provenance(q, "tsv")
    assert len(lines) == 2 and all("\n" not in line and "\t" not in line for line in lines)
