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


def test_committed_openapi_schema_matches_the_app() -> None:
    """Run ``python3 tools/export_openapi.py`` after changing a route."""
    import json
    from pathlib import Path

    from main import app
    committed = json.loads((Path(__file__).resolve().parents[1] / "openapi.json").read_text())
    assert committed == json.loads(json.dumps(app.openapi())), "stale api/openapi.json"


def test_annotation_genes_cursor_walks_the_same_pages_as_offset(client: TestClient) -> None:
    first = client.get("/annotations/ann-human/genes", params={"limit": 1}).json()
    assert first["next_cursor"] and first["offset"] == 0
    second = client.get("/annotations/ann-human/genes",
                        params={"limit": 1, "cursor": first["next_cursor"]}).json()
    by_offset = client.get("/annotations/ann-human/genes", params={"limit": 1, "offset": 1}).json()
    assert second["results"] == by_offset["results"] and second["offset"] == 1
    assert second["next_cursor"] is None                 # 2 genes in total: this was the last page
    whole = client.get("/annotations/ann-human/genes", params={"limit": 200}).json()
    assert whole["next_cursor"] is None


def test_annotation_cursor_errors_are_422(client: TestClient) -> None:
    both = client.get("/annotations/ann-human/genes",
                      params={"limit": 1, "offset": 1, "cursor": "eyJvIjoxfQ"})
    assert both.status_code == 422 and both.json()["detail"]["code"] == "invalid_combination"
    bad = client.get("/annotations/ann-human/genes", params={"cursor": "###"})
    assert bad.status_code == 422 and bad.json()["detail"]["code"] == "invalid_cursor"


def test_namespaces_list_has_no_next_page(client: TestClient) -> None:
    body = client.get("/annotations/ann-human/namespaces").json()
    assert body["next_cursor"] is None and body["limit"] == body["total"]


def test_accessions_prefix_is_the_same_filter_as_q(client: TestClient) -> None:
    path = "/annotations/ann-human/namespaces/refseq_transcript/accessions"
    a = client.get(path, params={"q": "XM_"}).json()
    b = client.get(path, params={"prefix": "XM_"}).json()
    assert a == b and a["total"] == 1 and a["results"][0]["accession"] == "XM_015139305.2"
    same = client.get(path, params={"q": "XM_", "prefix": "XM_"})
    assert same.status_code == 200
    clash = client.get(path, params={"q": "XM_0", "prefix": "XM_1"})
    assert clash.status_code == 422 and clash.json()["detail"]["code"] == "invalid_combination"


def test_namespaces_have_kind_and_slug(client: TestClient) -> None:
    rows = {r["namespace"]: r for r in client.get("/annotations/ann-human/namespaces").json()["results"]}
    assert rows["symbol"]["kind"] == "name"
    assert rows["refseq_transcript"]["kind"] == "database"
    assert rows["uniprot/swiss-prot"]["slug"] == "uniprot-swiss-prot"
    assert rows["symbol"]["slug"] == "symbol"
    r = client.get("/annotations/ann-human/namespaces/uniprot-swiss-prot/accessions")
    assert r.status_code == 200 and r.json()["total"] == 1


def test_genes_by_locus_ref(client: TestClient) -> None:
    nested = client.get("/annotations/ann-human/genes/10").json()
    assert client.get("/genes/ann-human:10").json() == nested
    assert client.get("/genes/ann-human:999").status_code == 404
    x1 = client.get("/annotations/ann-human/genes/10/xrefs").json()
    assert client.get("/genes/ann-human:10/xrefs").json() == x1
    for bad in ("ann-human", "ann-human:x", ":10", "ann-human:-1"):
        r = client.get(f"/genes/{bad}")
        assert r.status_code == 422, bad
        assert r.json()["detail"]["code"] == "invalid_locus_ref"


def test_locus_ref_with_unicode_digits_or_a_huge_id_is_422_not_500(client) -> None:
    for ref in ("ann-human:\u00b2", "ann-human:99999999999999999999999", "ann-human:-1", "nocolon"):
        r = client.get(f"/genes/{ref}")
        assert r.status_code == 422 and r.json()["detail"]["code"] == "invalid_locus_ref", ref


def test_a_name_shared_by_more_genes_than_sqlite_can_bind_still_pages(client, corpus_dir) -> None:
    import sqlite3

    n = 260_000                          # above SQLITE_MAX_VARIABLE_NUMBER (250000 in the image)
    conn = sqlite3.connect(corpus_dir / "per_gff" / "ann-human" / "genes.sqlite")
    conn.executemany(
        "INSERT INTO gene (local_id, source_gene_id, feature_type, seqid, start, end, strand) "
        "VALUES (?, ?, 'gene', 'chr1', ?, ?, 1)", [(1000 + i, f"g{i}", i, i + 5) for i in range(n)])
    conn.executemany("INSERT INTO gene_xref VALUES (?, 'alias', 'bigname', NULL, 'attr', 0)",
                     [(1000 + i,) for i in range(n)])
    conn.commit()
    conn.close()
    first = client.get("/annotations/ann-human/genes", params={"names": "bigname", "limit": 200})
    assert first.status_code == 200 and first.json()["total"] == n
    assert [g["local_id"] for g in first.json()["results"]][:2] == [1000, 1001]
    assert first.json()["next_cursor"]
    last = client.get("/annotations/ann-human/genes",
                      params={"names": "bigname", "limit": 200, "offset": 10_000}).json()
    assert last["total"] == n and len(last["results"]) == 200
    assert last["results"][0]["matches"][0]["name"] == "bigname"
