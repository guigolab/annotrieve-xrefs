"""API integration tests for the reshaped v1 surface."""
from __future__ import annotations

from fastapi.testclient import TestClient


def test_health(client: TestClient) -> None:
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_removed_top_level_routes(client: TestClient) -> None:
    assert client.get("/namespaces").status_code == 404
    assert client.get("/matches", params={"q": "tp53"}).status_code == 404
    assert client.get("/annotations").status_code == 404
    assert client.post(
        "/annotations", json={"curies": ["symbol:tp53"]}
    ).status_code == 404
    assert client.get("/annotations/index").status_code == 404


def test_gene_detail(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10")
    assert r.status_code == 200
    body = r.json()
    assert body["local_id"] == 10
    assert body["primary_name"] == "TP53"
    assert body["strand"] == "-"


def test_gene_detail_missing(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/999")
    assert r.status_code == 404
    r2 = client.get("/annotations/missing/genes/10")
    assert r2.status_code == 404


def test_annotation_namespaces(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/namespaces")
    assert r.status_code == 200
    body = r.json()
    ns = {row["namespace"] for row in body["results"]}
    assert "symbol" in ns
    assert "ensembl_transcript" in ns
    assert "offset" in body


def test_annotation_genes_browse(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 2
    assert body["offset"] == 0
    assert {row["local_id"] for row in body["results"]} == {10, 20}
    assert "next" not in body


def test_annotation_genes_offset(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"limit": 1, "offset": 0},
    )
    assert r.status_code == 200
    assert r.json()["total"] == 2
    assert len(r.json()["results"]) == 1
    first = r.json()["results"][0]["local_id"]
    r2 = client.get(
        "/annotations/ann-human/genes",
        params={"limit": 1, "offset": 1},
    )
    assert r2.json()["results"][0]["local_id"] != first


def test_annotation_genes_q(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"q": "tp"},
    )
    assert r.status_code == 200
    assert [row["primary_name"] for row in r.json()["results"]] == ["TP53"]


def test_annotation_genes_curies(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53"},
    )
    assert r.status_code == 200
    assert r.json()["total"] == 1
    assert r.json()["results"][0]["local_id"] == 10


def test_annotation_genes_curies_match_all(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,go:GO:0008150", "match": "all"},
    )
    assert r.status_code == 200
    assert r.json()["total"] == 1
    empty = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,symbol:brca1", "match": "all"},
    )
    assert empty.json()["total"] == 0


def test_gene_xrefs_list(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10/xrefs")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] >= 3
    assert body["offset"] == 0
    ns = {row["namespace"] for row in body["results"]}
    assert "symbol" in ns
    assert "go" in ns


def test_gene_xrefs_offset(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes/10/xrefs",
        params={"limit": 1, "offset": 0},
    )
    assert r.status_code == 200
    assert len(r.json()["results"]) == 1
    assert r.json()["total"] >= 3


def test_gene_xrefs_missing(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/999/xrefs")
    assert r.status_code == 404
