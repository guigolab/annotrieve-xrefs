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
    assert "annotation_id" not in body


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
    assert "refseq_transcript" in ns
    assert "offset" in body
    assert "annotation_id" not in body


def test_namespace_accessions_browse(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/namespaces/symbol/accessions")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 2
    assert "annotation_id" not in body
    by_acc = {row["accession"]: row for row in body["results"]}
    assert by_acc["tp53"]["display"] == "TP53"
    assert by_acc["tp53"]["gene_count"] == 1
    assert by_acc["brca1"]["display"] == "BRCA1"


def test_namespace_accessions_prefix(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/namespaces/refseq_transcript/accessions",
        params={"q": "XM_"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    row = body["results"][0]
    assert row["accession"] == "XM_015139305.2"
    assert row["display"] == "XM_015139305.2"
    assert row["gene_count"] == 1


def test_namespace_accessions_missing(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/namespaces/not_a_namespace/accessions"
    )
    assert r.status_code == 404
    assert r.json()["detail"]["code"] == "namespace_not_found"


def test_namespace_accessions_q_too_short(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/namespaces/symbol/accessions",
        params={"q": "ab"},
    )
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "prefix_too_short"


def test_namespace_accessions_symbol_q_casefold(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/namespaces/symbol/accessions",
        params={"q": "TP5"},
    )
    assert r.status_code == 200
    assert r.json()["total"] == 1
    assert r.json()["results"][0]["accession"] == "tp53"


def test_namespace_accessions_slash_and_short_prefix(
    client: TestClient,
) -> None:
    slash = client.get(
        "/annotations/ann-human/namespaces/uniprot/swiss-prot/accessions"
    )
    assert slash.status_code == 200
    assert slash.json()["total"] == 1
    assert slash.json()["results"][0]["accession"] == "P04637"
    assert slash.json()["results"][0]["display"] == "P53_HUMAN"

    encoded = client.get(
        "/annotations/ann-human/namespaces/uniprot%2Fswiss-prot/accessions"
    )
    assert encoded.status_code == 200
    assert encoded.json()["total"] == 1

    short = client.get(
        "/annotations/ann-human/namespaces/swiss-prot/accessions"
    )
    assert short.status_code == 200
    assert short.json()["total"] == 1
    assert short.json()["results"][0]["accession"] == "P04637"


def test_annotation_genes_browse(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 2
    assert body["offset"] == 0
    assert {row["local_id"] for row in body["results"]} == {10, 20}
    assert "next" not in body
    assert "annotation_id" not in body
    assert "matches" not in body["results"][0]
    assert "annotation_id" not in body["results"][0]


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


def test_annotation_genes_names_exact(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"names": "TP53"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert "annotation_id" not in body
    gene = body["results"][0]
    assert gene["local_id"] == 10
    assert "annotation_id" not in gene
    assert gene["matches"] == [
        {"input": "TP53", "name": "tp53", "ns": "alias"},
        {"input": "TP53", "name": "tp53", "ns": "symbol"},
    ]


def test_annotation_genes_names_rejects_wildcard(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"names": "tp5*"},
    )
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "wildcard_not_supported"


def test_annotation_genes_names_match_any_all(client: TestClient) -> None:
    any_r = client.get(
        "/annotations/ann-human/genes",
        params={"names": "tp53,brca1", "match": "any"},
    )
    assert any_r.status_code == 200
    assert any_r.json()["total"] == 2
    by_id = {g["local_id"]: g for g in any_r.json()["results"]}
    assert by_id[10]["matches"] == [
        {"input": "tp53", "name": "tp53", "ns": "alias"},
        {"input": "tp53", "name": "tp53", "ns": "symbol"},
    ]
    assert by_id[20]["matches"] == [
        {"input": "brca1", "name": "brca1", "ns": "symbol"},
    ]
    all_r = client.get(
        "/annotations/ann-human/genes",
        params={"names": "tp53,brca1", "match": "all"},
    )
    assert all_r.status_code == 200
    assert all_r.json()["total"] == 0


def test_annotation_genes_names_validation(client: TestClient) -> None:
    empty = client.get(
        "/annotations/ann-human/genes",
        params={"names": " , "},
    )
    assert empty.status_code == 422
    assert empty.json()["detail"]["code"] == "names_required"

    short = client.get(
        "/annotations/ann-human/genes",
        params={"names": "ab*"},
    )
    assert short.status_code == 422
    assert short.json()["detail"]["code"] == "wildcard_not_supported"

    bad = client.get(
        "/annotations/ann-human/genes",
        params={"names": "tp53", "match": "xor"},
    )
    assert bad.status_code == 422
    assert bad.json()["detail"]["code"] == "invalid_match"


def test_annotation_genes_curies_exact(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "refseq_transcript:XM_015139305.2"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    gene = body["results"][0]
    assert gene["local_id"] == 10
    assert gene["matches"] == [
        {
            "input": "refseq_transcript:XM_015139305.2",
            "ns": "refseq_transcript",
            "accession": "XM_015139305.2",
        }
    ]


def test_annotation_genes_curies_go_and_swiss_prot(client: TestClient) -> None:
    go = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "go:8150"},
    )
    assert go.status_code == 200
    assert go.json()["total"] == 1
    assert go.json()["results"][0]["matches"] == [
        {"input": "go:8150", "ns": "go", "accession": "GO:0008150"}
    ]

    sp = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "swiss-prot:P04637"},
    )
    assert sp.status_code == 200
    assert sp.json()["total"] == 1
    assert sp.json()["results"][0]["local_id"] == 10
    assert sp.json()["results"][0]["matches"] == [
        {
            "input": "swiss-prot:P04637",
            "ns": "uniprot/swiss-prot",
            "accession": "P04637",
        }
    ]


def test_annotation_genes_curies_match_any_all(client: TestClient) -> None:
    any_r = client.get(
        "/annotations/ann-human/genes",
        params={
            "curies": "symbol:tp53,symbol:brca1",
            "match": "any",
        },
    )
    assert any_r.status_code == 200
    assert any_r.json()["total"] == 2
    all_r = client.get(
        "/annotations/ann-human/genes",
        params={
            "curies": "symbol:tp53,symbol:brca1",
            "match": "all",
        },
    )
    assert all_r.status_code == 200
    assert all_r.json()["total"] == 0


def test_annotation_genes_curies_miss(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "refseq_transcript:XM_000000000.1"},
    )
    assert r.status_code == 200
    assert r.json()["total"] == 0


def test_annotation_genes_names_and_curies_conflict(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"names": "tp53", "curies": "symbol:tp53"},
    )
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "names_and_curies"


def test_annotation_genes_curies_malformed(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "tp53"},
    )
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "missing_colon"


def test_gene_xrefs_list(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10/xrefs")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] >= 3
    assert body["offset"] == 0
    assert "annotation_id" not in body
    assert "local_id" not in body
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
