"""HTTP-level checks for annotrieve-xrefs v1 routes."""
from __future__ import annotations

from fastapi.testclient import TestClient


def test_health(client: TestClient) -> None:
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_namespaces(client: TestClient) -> None:
    r = client.get("/namespaces")
    assert r.status_code == 200
    body = r.json()
    names = {row["namespace"] for row in body["results"]}
    assert names == {"alias", "go", "interpro", "symbol"}
    assert body["next"] is None


def test_namespace_item(client: TestClient) -> None:
    r = client.get("/namespaces/symbol")
    assert r.status_code == 200
    assert r.json() == {"namespace": "symbol", "accession_count": 1}
    missing = client.get("/namespaces/nope")
    assert missing.status_code == 404


def test_accessions_sorted_by_n_loci(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/accessions")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["results"][0]["accession"] == "tp53"
    assert body["results"][0]["n_annotations"] == 2
    assert body["results"][0]["n_loci"] == 2


def test_genes_any_match_merges_matched(client: TestClient) -> None:
    r = client.get(
        "/genes",
        params=[("id", "symbol:TP53"), ("id", "alias:tp53")],
    )
    assert r.status_code == 200
    body = r.json()
    assert "total" not in body
    assert len(body["results"]) == 2
    human = next(x for x in body["results"] if x["annotation_id"] == "ann-human")
    assert set(human["matched"]) == {"symbol:tp53", "alias:tp53"}
    assert human["local_id"] == 10
    mouse = next(x for x in body["results"] if x["annotation_id"] == "ann-mouse")
    assert mouse["matched"] == ["symbol:tp53"]


def test_genes_taxid_filters_lineage(client: TestClient) -> None:
    # Homo sapiens only
    r = client.get(
        "/genes",
        params={"id": "symbol:tp53", "taxid": 9606},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["taxid"] == 9606
    assert len(body["results"]) == 1
    assert body["results"][0]["annotation_id"] == "ann-human"

    # Mammalia keeps both
    r2 = client.get(
        "/genes",
        params={"id": "symbol:tp53", "taxid": 40674},
    )
    assert len(r2.json()["results"]) == 2


def test_genes_unknown_prefix_error_entry(client: TestClient) -> None:
    r = client.get(
        "/genes",
        params=[("id", "symbol:tp53"), ("id", "notans:foo")],
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    assert body["errors"]
    assert "unknown CURIE prefix" in body["errors"][0]["error"]


def test_genes_go_padding(client: TestClient) -> None:
    r = client.get("/genes", params={"id": "GO:8150"})
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    for row in body["results"]:
        assert row["matched"] == ["go:GO:0008150"]


def test_genes_interpro_case(client: TestClient) -> None:
    r = client.get("/genes", params={"id": "interpro:IPR002117"})
    assert r.status_code == 200
    assert len(r.json()["results"]) == 1


def test_genes_requires_id(client: TestClient) -> None:
    r = client.get("/genes")
    assert r.status_code == 422


def test_gene_detail(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10")
    assert r.status_code == 200
    body = r.json()
    assert body["annotation_id"] == "ann-human"
    assert body["local_id"] == 10
    assert body["prose"] == "tumor protein p53"
    ns = {x["namespace"] for x in body["xrefs"]}
    assert "symbol" in ns
    assert "go" in ns


def test_gene_detail_missing(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/999")
    assert r.status_code == 404
    r2 = client.get("/annotations/missing/genes/10")
    assert r2.status_code == 404


def test_genes_rejects_nonpositive_taxid(client: TestClient) -> None:
    r = client.get("/genes", params={"id": "symbol:tp53", "taxid": 0})
    assert r.status_code == 422


def test_genes_too_many_ids(client: TestClient) -> None:
    params = [("id", f"symbol:g{i}") for i in range(33)]
    r = client.get("/genes", params=params)
    assert r.status_code == 400
    assert "at most" in r.json()["detail"]


def test_gene_detail_rejects_negative_local_id(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/-1")
    assert r.status_code == 422


def test_oversized_cursor_rejected(client: TestClient) -> None:
    huge = "A" * 5000
    r = client.get("/genes", params={"id": "symbol:tp53", "cursor": huge})
    assert r.status_code == 400

    first = client.get(
        "/genes",
        params={"id": "symbol:tp53", "limit": 1},
    ).json()
    assert len(first["results"]) == 1
    assert first["next"]
    second = client.get(
        "/genes",
        params={
            "id": "symbol:tp53",
            "limit": 1,
            "cursor": first["next"],
        },
    ).json()
    assert len(second["results"]) == 1
    assert first["results"][0]["annotation_id"] != second["results"][0]["annotation_id"]
    assert second["previous"]
