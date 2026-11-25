"""/distributions over a small sidecar-built index with known gene lengths and two clades."""
from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from gene_corpus.serve_build_cli import main as serve_build_main
from gene_corpus.tests.test_catalogue_backfill import tsv
from gene_corpus.tests.test_sidecar_manifest import make_shard

MAMMALS, BIRDS = 40674, 8782


def gene(local_id: int, start: int, end: int, name: str | None = None):
    return (local_id, f"gene-{local_id}", "gene", "chr1", start, end, 1, None, name, None, None, 0)


def sym(local_id: int, name: str):
    return (local_id, "symbol", name, None, "gff", 0)


def alias(local_id: int, name: str):
    return (local_id, "alias", name, None, "gff", 0)


# annotation, taxid, source, lineage, genes, xrefs      (lengths are end - start + 1)
ANNOTATIONS = [
    ("m0", 9606, "RefSeq", f"1;{MAMMALS};9606", [gene(1, 1, 1000)], [sym(1, "tp53")]),
    # two symbol copies: 2000 and 4000 -> mean 3000, copies 2
    ("m1", 10090, "RefSeq", f"1;{MAMMALS};10090", [gene(1, 1, 2000), gene(2, 5000, 8999)],
     [sym(1, "tp53"), sym(2, "tp53")]),
    # only an alias carries the name here: 500 bp
    ("m2", 9913, "Ensembl", f"1;{MAMMALS};9913", [gene(1, 1, 500)], [alias(1, "tp53"), sym(1, "other")]),
    ("b0", 9031, "RefSeq", f"1;{BIRDS};9031", [gene(1, 1, 10000)], [sym(1, "tp53")]),
    ("b1", 9103, "RefSeq", f"1;{BIRDS};9103", [gene(1, 1, 300)], [sym(1, "other")]),
]


@pytest.fixture()
def dist_client(tmp_path: Path):
    import settings as settings_mod
    from helpers import presence as presence_mod
    from main import create_app

    work = tmp_path / "corpus"
    for aid, _t, source, _l, genes, xrefs in ANNOTATIONS:
        make_shard(work / "per_gff", aid, source=source, genes=genes, xrefs=xrefs)
    rep = tmp_path / "report.tsv"
    rep.write_text("\n".join(tsv(*[(aid, "GCF_1", "n", "o", taxid, source, "p", "2024-01-01", lin)
                                   for aid, taxid, source, lin, _g, _x in ANNOTATIONS])) + "\n")
    assert serve_build_main(["sidecars", "--work-dir", str(work), "--workers", "1"]) == 0
    assert serve_build_main(["catalogue", "--work-dir", str(work), "--tsv-file", str(rep)]) == 0
    assert serve_build_main(["build", "--work-dir", str(work)]) == 0
    assert serve_build_main(["loci", "--work-dir", str(work)]) == 0
    old = settings_mod.settings.PRESENCE_DIR
    settings_mod.settings.PRESENCE_DIR = str(work / "serve")
    presence_mod.close_presence()
    with TestClient(create_app()) as client:
        yield client
    presence_mod.close_presence()
    settings_mod.settings.PRESENCE_DIR = old


COLUMNS = [{"id": "mammals", "taxids": [MAMMALS]}, {"id": "birds", "taxids": [BIRDS]}]


def cells(body: dict) -> dict[tuple[str, str], dict]:
    names = [n["name"] for n in body["names"]]
    cols = [c["id"] for c in body["columns"]]
    return {(names[r[0]], cols[r[1]]): dict(zip(body["fields"], r)) for r in body["rows"]}


def test_summary_per_clade_with_exact_stats(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS}).json()
    assert body["fields"] == ["name", "column", "n_with_gene", "n_multi_copy", "mean", "min", "p25",
                              "median", "p75", "max", "hist", "p5", "p95"]
    assert [c["n_annotations"] for c in body["columns"]] == [3, 2]
    got = cells(body)
    m, b = got[("tp53", "mammals")], got[("tp53", "birds")]
    # mammals: 1000 (m0), 3000 (mean of 2000 and 4000, m1), 500 (alias-only gene, m2)
    assert (m["n_with_gene"], m["n_multi_copy"]) == (3, 1)
    assert (m["min"], m["median"], m["max"]) == (500.0, 1000.0, 3000.0)
    assert m["mean"] == 1500.0 and m["p25"] == 750.0 and m["p75"] == 2000.0
    # whiskers: type-7 quantiles of [500, 1000, 3000]; one value gives itself
    assert m["p5"] == 550.0 and m["p95"] == 2800.0
    assert b["p5"] == b["p95"] == 10000.0
    assert (b["n_with_gene"], b["n_multi_copy"], b["median"]) == (1, 0, 10000.0)
    assert body["bins"] == {"scale": "log10", "min_bp": 100, "max_bp": 10_000_000, "n": 40}
    assert sum(m["hist"]) == 3 and len(m["hist"]) == 40 and sum(b["hist"]) == 1
    assert body["total"] == len(body["rows"]) == 2


def test_alias_off_drops_alias_only_annotations(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS,
                                                    "alias": "off"}).json()
    m = cells(body)[("tp53", "mammals")]
    assert m["n_with_gene"] == 2 and m["min"] == 1000.0


def test_cells_without_the_gene_are_left_out_and_names_expand(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["other", "nosuchgene"],
                                                    "columns": COLUMNS}).json()
    got = cells(body)
    assert set(got) == {("other", "mammals"), ("other", "birds")}
    assert got[("other", "birds")]["median"] == 300.0
    assert [t["found"] for t in body["terms"]] == [True, False]


def test_overlapping_columns_and_default_column(dist_client) -> None:
    both = dist_client.post("/distributions", json={"names": ["tp53"], "columns": [
        {"id": "everything", "taxids": [1]}, {"id": "mammals", "taxids": [MAMMALS]}]}).json()
    got = cells(both)
    assert got[("tp53", "everything")]["n_with_gene"] == 4
    assert got[("tp53", "mammals")]["n_with_gene"] == 3
    default = dist_client.post("/distributions", json={"names": ["tp53"]}).json()
    assert [c["id"] for c in default["columns"]] == ["all"]
    assert cells(default)[("tp53", "all")]["n_with_gene"] == 4


def test_filters_scope_the_denominator_and_the_cells(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS,
                                                    "sources": ["Ensembl"]}).json()
    assert [c["n_annotations"] for c in body["columns"]] == [1, 0]
    assert set(cells(body)) == {("tp53", "mammals")}
    ex = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS,
                                                  "exclude_taxids": [BIRDS]}).json()
    assert [c["n_annotations"] for c in ex["columns"]] == [3, 0]


def test_prefix_term_gives_one_row_group_per_expanded_name(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["tp5*"], "columns": COLUMNS}).json()
    assert [n["name"] for n in body["names"]] == ["tp53"]
    assert body["terms"][0]["kind"] == "prefix"


def test_dry_run_counts_only(dist_client) -> None:
    body = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS,
                                                    "dry_run": True}).json()
    assert "rows" not in body and body["n_annotations"] == 4
    assert body["columns"][0]["n_annotations"] == 3


def test_errors(dist_client) -> None:
    post = dist_client.post
    assert post("/distributions", json={"names": ["tp53"], "columns": []}).status_code == 422
    assert post("/distributions", json={"names": ["tp53"], "columns": [
        {"id": "a", "taxids": [1]}, {"id": "a", "taxids": [2]}]}).status_code == 422
    assert post("/distributions", json={"names": ["tp53"], "columns": [
        {"id": str(i), "taxids": [1]} for i in range(51)]}).status_code == 422
    assert post("/distributions", json={"names": ["tp53"], "alias": "x"}).status_code == 422
    assert post("/distributions", json={"columns": COLUMNS}).status_code == 422
    assert post("/distributions", json={"names": ["tp53"], "bogus": 1}).status_code == 422


def test_too_large_request_is_refused(dist_client) -> None:
    import settings as settings_mod

    old = settings_mod.settings.DISTRIBUTIONS_MAX_LOCI
    settings_mod.settings.DISTRIBUTIONS_MAX_LOCI = 2
    try:
        r = dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS})
        assert r.status_code == 413 and r.json()["detail"]["code"] == "distributions_too_large"
        assert dist_client.post("/distributions", json={"names": ["tp53"], "dry_run": True}).status_code == 200
    finally:
        settings_mod.settings.DISTRIBUTIONS_MAX_LOCI = old


def test_points_of_one_clade(dist_client) -> None:
    body = dist_client.get(f"/distributions/{MAMMALS}/points", params={"name": "tp53"}).json()
    assert body["fields"] == ["annotation_id", "taxid", "length", "copies", "via"]
    assert body["taxid"] == MAMMALS and body["name"] == "tp53" and body["n_annotations"] == 3
    # (the shard fixture stores taxid 9606 for every annotation)
    assert body["rows"] == [["m0", 9606, 1000.0, 1, 0], ["m1", 9606, 3000.0, 2, 0],
                            ["m2", 9606, 500.0, 1, 1]]
    assert body["total"] == 3
    off = dist_client.get(f"/distributions/{MAMMALS}/points", params={"name": "tp53", "alias": "off"}).json()
    assert [r[0] for r in off["rows"]] == ["m0", "m1"]
    birds = dist_client.get(f"/distributions/{BIRDS}/points", params={"name": "tp53"}).json()
    assert birds["rows"] == [["b0", 9606, 10000.0, 1, 0]]


def test_points_agree_with_the_summary(dist_client) -> None:
    summary = cells(dist_client.post("/distributions", json={"names": ["tp53"], "columns": COLUMNS}).json())
    for taxid, col in ((MAMMALS, "mammals"), (BIRDS, "birds")):
        pts = dist_client.get(f"/distributions/{taxid}/points", params={"name": "tp53"}).json()
        lengths = sorted(r[2] for r in pts["rows"])
        assert len(lengths) == summary[("tp53", col)]["n_with_gene"]
        assert lengths[0] == summary[("tp53", col)]["min"] and lengths[-1] == summary[("tp53", col)]["max"]


def test_points_errors_and_empty(dist_client) -> None:
    get = dist_client.get
    assert get(f"/distributions/{MAMMALS}/points", params={"name": "tp5*"}).status_code == 422
    assert get(f"/distributions/{MAMMALS}/points", params={"name": "tp53,other"}).status_code == 422
    assert get(f"/distributions/{MAMMALS}/points").status_code == 422
    assert get(f"/distributions/{MAMMALS}/points", params={"name": "tp53", "x": 1}).status_code == 422
    assert get("/distributions/0/points", params={"name": "tp53"}).status_code == 422
    none = get(f"/distributions/{MAMMALS}/points", params={"name": "nosuchgene"}).json()
    assert none["rows"] == [] and none["terms"][0]["found"] is False
    empty = get("/distributions/999999/points", params={"name": "tp53"}).json()
    assert empty["rows"] == [] and empty["n_annotations"] == 0


def test_points_are_cacheable_by_data_version(dist_client) -> None:
    r = dist_client.get(f"/distributions/{BIRDS}/points", params={"name": "tp53"})
    etag = r.headers.get("etag")
    assert etag and "max-age" in r.headers.get("cache-control", "")
    again = dist_client.get(f"/distributions/{BIRDS}/points", params={"name": "tp53"},
                            headers={"If-None-Match": etag})
    assert again.status_code == 304


def test_bin_and_quantile_helpers() -> None:
    from services.distributions import BIN_COUNT, bin_of, quantile

    assert bin_of(1) == 0 and bin_of(100) == 0 and bin_of(10_000_000) == BIN_COUNT - 1
    assert bin_of(1e9) == BIN_COUNT - 1
    assert bin_of(1_000) == 8 and bin_of(100_000) == 24         # 0.125 decade per bucket
    assert quantile([5.0], 0.5) == 5.0
    assert quantile([1.0, 2.0, 3.0, 4.0], 0.5) == 2.5
    assert quantile([1.0, 2.0, 3.0, 4.0], 0.25) == 1.75
