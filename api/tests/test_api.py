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
    assert missing.json()["detail"] == {
        "code": "namespace_not_found",
        "message": "Unknown namespace: nope",
    }


def test_accessions_sorted_by_n_genes(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/accessions")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["results"][0]["accession"] == "tp53"
    assert body["results"][0]["n_annotations"] == 2
    assert body["results"][0]["n_genes"] == 2
    assert "n_loci" not in body["results"][0]


def test_accessions_filter_casefold_and_missing(client: TestClient) -> None:
    r = client.get(
        "/namespaces/symbol/accessions",
        params={"accessions": "tp53,TP53,missing"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert len(body["results"]) == 1
    assert body["results"][0]["accession"] == "tp53"


def test_accessions_filter_go_canonical(client: TestClient) -> None:
    r = client.get(
        "/namespaces/go/accessions",
        params={"accessions": "8150,GO:8150"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert len(body["results"]) == 1
    assert body["results"][0]["accession"] == "GO:0008150"


def test_accessions_filter_too_many(client: TestClient) -> None:
    vals = ",".join(f"a{i}" for i in range(21))
    r = client.get(
        "/namespaces/symbol/accessions",
        params={"accessions": vals},
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "too_many_accessions"


def test_accession_detail(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/accessions/tp53")
    assert r.status_code == 200
    assert r.json() == {
        "namespace": "symbol",
        "accession": "tp53",
        "n_annotations": 2,
        "n_genes": 2,
    }
    folded = client.get("/namespaces/symbol/accessions/TP53")
    assert folded.status_code == 200
    assert folded.json()["accession"] == "tp53"
    missing = client.get("/namespaces/symbol/accessions/nope")
    assert missing.status_code == 404
    assert missing.json()["detail"]["code"] == "accession_not_found"
    bad_ns = client.get("/namespaces/nope/accessions/tp53")
    assert bad_ns.status_code == 404
    assert bad_ns.json()["detail"]["code"] == "namespace_not_found"


def test_accession_detail_go(client: TestClient) -> None:
    for path in (
        "/namespaces/go/accessions/8150",
        "/namespaces/go/accessions/GO:0008150",
        "/namespaces/go/accessions/GO:8150",
    ):
        r = client.get(path)
        assert r.status_code == 200, path
        body = r.json()
        assert body["namespace"] == "go"
        assert body["accession"] == "GO:0008150"
        assert body["n_annotations"] == 2
        assert body["n_genes"] == 2


def test_hit_annotations_collapse(client: TestClient) -> None:
    r = client.get(
        "/hits/annotations",
        params={"curies": "symbol:tp53"},
    )
    assert r.status_code == 200
    body = r.json()
    assert "total" not in body
    assert len(body["results"]) == 2
    by_id = {row["annotation_id"]: row for row in body["results"]}
    assert by_id["ann-human"]["n_genes"] == 1
    assert by_id["ann-human"]["matched"] == ["symbol:tp53"]
    assert by_id["ann-mouse"]["n_genes"] == 1
    assert "local_id" not in by_id["ann-human"]
    assert "example" not in by_id["ann-human"]


def test_hit_annotations_taxid(client: TestClient) -> None:
    r = client.get(
        "/hits/annotations",
        params={"curies": "symbol:tp53", "taxid": 9606},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["taxid"] == 9606
    assert len(body["results"]) == 1
    assert body["results"][0]["annotation_id"] == "ann-human"


def test_hit_annotations_multi_curie_no_double_count(client: TestClient) -> None:
    r = client.get(
        "/hits/annotations",
        params={"curies": "symbol:tp53,alias:tp53"},
    )
    assert r.status_code == 200
    body = r.json()
    human = next(x for x in body["results"] if x["annotation_id"] == "ann-human")
    assert human["n_genes"] == 1
    assert set(human["matched"]) == {"symbol:tp53", "alias:tp53"}


def test_hit_annotations_post(client: TestClient) -> None:
    r = client.post(
        "/hits/annotations",
        json={"curies": ["symbol:tp53"], "taxid": 9606},
    )
    assert r.status_code == 200
    assert len(r.json()["results"]) == 1
    assert r.json()["results"][0]["annotation_id"] == "ann-human"


def test_hits_any_match_merges_matched(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:TP53,alias:tp53"},
    )
    assert r.status_code == 200
    body = r.json()
    assert "total" not in body
    assert len(body["results"]) == 2
    human = next(x for x in body["results"] if x["annotation_id"] == "ann-human")
    assert set(human["matched"]) == {"symbol:tp53", "alias:tp53"}
    assert human["local_id"] == 10
    assert human["strand"] == "-"
    mouse = next(x for x in body["results"] if x["annotation_id"] == "ann-mouse")
    assert mouse["matched"] == ["symbol:tp53"]
    assert mouse["strand"] == "+"


def test_hits_taxid_filters_lineage(client: TestClient) -> None:
    # Homo sapiens only
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "taxid": 9606},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["taxid"] == 9606
    assert len(body["results"]) == 1
    assert body["results"][0]["annotation_id"] == "ann-human"

    # Mammalia keeps both
    r2 = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "taxid": 40674},
    )
    assert len(r2.json()["results"]) == 2


def test_hits_unknown_prefix_error_entry(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,notans:foo"},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    assert body["errors"] == [
        {
            "curie": "notans:foo",
            "code": "unknown_prefix",
            "message": "unknown CURIE prefix: notans",
        }
    ]


def test_hits_go_padding(client: TestClient) -> None:
    r = client.get("/hits", params={"curies": "GO:8150"})
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    for row in body["results"]:
        assert row["matched"] == ["go:GO:0008150"]


def test_hits_interpro_case(client: TestClient) -> None:
    r = client.get("/hits", params={"curies": "interpro:IPR002117"})
    assert r.status_code == 200
    assert len(r.json()["results"]) == 1


def test_hits_requires_curies(client: TestClient) -> None:
    r = client.get("/hits")
    assert r.status_code == 422
    detail = r.json()["detail"]
    assert detail["code"] == "validation_error"
    assert "curies" in detail["message"]
    assert isinstance(detail["errors"], list)
    assert detail["errors"]
    missing = client.post("/hits", json={})
    assert missing.status_code == 422
    assert missing.json()["detail"]["code"] == "validation_error"


def test_gene_detail(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10")
    assert r.status_code == 200
    body = r.json()
    assert body["annotation_id"] == "ann-human"
    assert body["local_id"] == 10
    assert body["strand"] == "-"
    assert body["prose"] == "tumor protein p53"
    assert body["primary_name"] == "TP53"
    assert "xrefs" not in body
    assert "matched" not in body


def test_gene_detail_missing(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/999")
    assert r.status_code == 404
    assert r.json()["detail"] == {
        "code": "gene_not_found",
        "message": "Gene local_id 999 not found in ann-human",
    }
    r2 = client.get("/annotations/missing/genes/10")
    assert r2.status_code == 404
    assert r2.json()["detail"]["code"] == "shard_not_found"


def test_hits_rejects_nonpositive_taxid(client: TestClient) -> None:
    r = client.get("/hits", params={"curies": "symbol:tp53", "taxid": 0})
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "validation_error"


def test_hits_too_many_curies(client: TestClient) -> None:
    curies = ",".join(f"symbol:g{i}" for i in range(21))
    r = client.get("/hits", params={"curies": curies})
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert detail["code"] == "too_many_curies"
    assert detail["max"] == 20
    assert detail["received"] == 21
    assert "at most" in detail["message"]


def test_hits_cursor_filter_mismatch_taxid(client: TestClient) -> None:
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "taxid": 9606, "limit": 1},
    )
    assert first.status_code == 200
    # Only one human hit → often no next; use mammalia page with 2 hits.
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "limit": 1},
    )
    assert first.status_code == 200
    token = first.json()["next"]
    assert token
    bad = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "taxid": 9606, "next": token},
    )
    assert bad.status_code == 400
    assert bad.json()["detail"]["code"] == "cursor_filter_mismatch"


def test_hits_cursor_same_filters_reordered_curies(client: TestClient) -> None:
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53,alias:tp53", "limit": 1},
    )
    assert first.status_code == 200
    token = first.json()["next"]
    assert token
    ok = client.get(
        "/hits",
        params={"curies": "alias:tp53,symbol:tp53", "limit": 1, "next": token},
    )
    assert ok.status_code == 200


def test_accessions_cursor_filter_mismatch(client: TestClient) -> None:
    from helpers.cursor import encode_accessions_cursor, filter_fingerprint

    filter_f = filter_fingerprint(
        {"sort": "n_genes", "sort_order": "desc", "accessions": ["tp53"]}
    )
    token = encode_accessions_cursor(
        sort="n_genes",
        sort_order="desc",
        sort_value=2,
        accession="tp53",
        filter_f=filter_f,
    )
    r = client.get(
        "/namespaces/symbol/accessions",
        params={"next": token},
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "cursor_filter_mismatch"

    ok = client.get(
        "/namespaces/symbol/accessions",
        params={"accessions": "TP53", "next": token},
    )
    assert ok.status_code == 200


def test_annotation_genes_cursor_filter_mismatch(client: TestClient) -> None:
    first = client.get(
        "/annotations/ann-human/genes",
        params={"q": "TP", "limit": 1},
    )
    assert first.status_code == 200
    token = first.json().get("next")
    if not token:
        # Single prefix match — craft via second page of unfiltered list.
        first = client.get(
            "/annotations/ann-human/genes",
            params={"limit": 1},
        )
        token = first.json()["next"]
        assert token
        bad = client.get(
            "/annotations/ann-human/genes",
            params={"q": "TP", "limit": 1, "next": token},
        )
    else:
        bad = client.get(
            "/annotations/ann-human/genes",
            params={"q": "symbol:tp53", "limit": 1, "next": token},
        )
    assert bad.status_code == 400
    assert bad.json()["detail"]["code"] == "cursor_filter_mismatch"


def test_hits_busy_returns_503(client: TestClient, monkeypatch) -> None:
    import helpers.concurrency as concurrency

    # Saturate the gate so the next acquire fails immediately.
    for _ in range(concurrency.HITS_MAX_INFLIGHT):
        assert concurrency._HITS_SEMAPHORE.acquire(blocking=False)

    try:
        r = client.get("/hits", params={"curies": "symbol:tp53"})
        assert r.status_code == 503
        assert r.headers.get("Retry-After") == "1"
        assert r.json()["detail"]["code"] == "hits_busy"
    finally:
        for _ in range(concurrency.HITS_MAX_INFLIGHT):
            concurrency._HITS_SEMAPHORE.release()


def test_gene_detail_rejects_negative_local_id(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/-1")
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "validation_error"


def test_oversized_page_token_rejected(client: TestClient) -> None:
    huge = "A" * 5000
    r = client.get("/hits", params={"curies": "symbol:tp53", "next": huge})
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "invalid_cursor"

    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "limit": 1},
    ).json()
    assert len(first["results"]) == 1
    assert first["next"]
    second = client.get(
        "/hits",
        params={
            "curies": "symbol:tp53",
            "limit": 1,
            "next": first["next"],
        },
    ).json()
    assert len(second["results"]) == 1
    assert first["results"][0]["annotation_id"] != second["results"][0]["annotation_id"]
    assert second["previous"]


def test_conflicting_pagination_params(client: TestClient) -> None:
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "limit": 1},
    ).json()
    assert first["next"]
    r = client.get(
        "/hits",
        params={
            "curies": "symbol:tp53",
            "limit": 1,
            "next": first["next"],
            "previous": first["next"],
        },
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "conflicting_pagination"


def test_hits_curies_csv_strips_empties(client: TestClient) -> None:
    r = client.get("/hits", params={"curies": " symbol:tp53 , "})
    assert r.status_code == 200
    assert len(r.json()["results"]) == 2


def test_hits_empty_curies(client: TestClient) -> None:
    r = client.get("/hits", params={"curies": " , "})
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "curies_required"
    posted = client.post("/hits", json={"curies": ["", "  "]})
    assert posted.status_code == 400
    assert posted.json()["detail"]["code"] == "curies_required"


def test_hits_post_matches_get(client: TestClient) -> None:
    r = client.post(
        "/hits",
        json={"curies": ["symbol:TP53", "alias:tp53"]},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    human = next(x for x in body["results"] if x["annotation_id"] == "ann-human")
    assert set(human["matched"]) == {"symbol:tp53", "alias:tp53"}


def test_hits_post_allows_more_than_get_limit(client: TestClient) -> None:
    ok = client.post(
        "/hits",
        json={"curies": [f"symbol:g{i}" for i in range(21)]},
    )
    assert ok.status_code == 200
    too_many = client.post(
        "/hits",
        json={"curies": [f"symbol:g{i}" for i in range(101)]},
    )
    assert too_many.status_code == 400
    detail = too_many.json()["detail"]
    assert detail["code"] == "too_many_curies"
    assert detail["max"] == 100
    assert detail["received"] == 101


def test_hits_match_all_same_gene(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,alias:tp53", "match": "all"},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 1
    row = body["results"][0]
    assert row["annotation_id"] == "ann-human"
    assert row["local_id"] == 10
    assert set(row["matched"]) == {"symbol:tp53", "alias:tp53"}


def test_hits_match_any_still_union(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,alias:tp53", "match": "any"},
    )
    assert r.status_code == 200
    assert len(r.json()["results"]) == 2


def test_hits_match_all_shared_go(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,GO:0008150", "match": "all"},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    by_id = {row["annotation_id"]: row for row in body["results"]}
    assert by_id["ann-human"]["local_id"] == 10
    assert by_id["ann-mouse"]["local_id"] == 3
    assert set(by_id["ann-human"]["matched"]) == {"go:GO:0008150", "symbol:tp53"}


def test_hit_annotations_match_all(client: TestClient) -> None:
    r = client.get(
        "/hits/annotations",
        params={"curies": "symbol:tp53,alias:tp53", "match": "all"},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 1
    assert body["results"][0]["annotation_id"] == "ann-human"
    assert body["results"][0]["n_genes"] == 1
    assert set(body["results"][0]["matched"]) == {"symbol:tp53", "alias:tp53"}


def test_hit_annotations_match_any_still_union(client: TestClient) -> None:
    r = client.get(
        "/hits/annotations",
        params={"curies": "symbol:tp53,alias:tp53", "match": "any"},
    )
    assert r.status_code == 200
    assert len(r.json()["results"]) == 2


def test_hits_match_all_missing_accession_empty(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,symbol:nosuchgene", "match": "all"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["results"] == []
    assert "errors" not in body or body.get("errors") is None


def test_hits_match_all_invalid_prefix_empty(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53,notans:foo", "match": "all"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["results"] == []
    assert body["errors"] == [
        {
            "curie": "notans:foo",
            "code": "unknown_prefix",
            "message": "unknown CURIE prefix: notans",
        }
    ]


def test_hits_invalid_match(client: TestClient) -> None:
    r = client.get(
        "/hits",
        params={"curies": "symbol:tp53", "match": "nope"},
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "invalid_match"


def test_hits_match_cursor_filter_mismatch(client: TestClient) -> None:
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53,GO:0008150", "match": "any", "limit": 1},
    )
    assert first.status_code == 200
    token = first.json()["next"]
    assert token
    bad = client.get(
        "/hits",
        params={
            "curies": "symbol:tp53,GO:0008150",
            "match": "all",
            "limit": 1,
            "next": token,
        },
    )
    assert bad.status_code == 400
    assert bad.json()["detail"]["code"] == "cursor_filter_mismatch"


def test_hits_match_all_pagination(client: TestClient) -> None:
    first = client.get(
        "/hits",
        params={"curies": "symbol:tp53,GO:0008150", "match": "all", "limit": 1},
    )
    assert first.status_code == 200
    body = first.json()
    assert len(body["results"]) == 1
    assert body["results"][0]["annotation_id"] == "ann-human"
    token = body["next"]
    assert token
    second = client.get(
        "/hits",
        params={
            "curies": "symbol:tp53,GO:0008150",
            "match": "all",
            "limit": 1,
            "next": token,
        },
    )
    assert second.status_code == 200
    page2 = second.json()
    assert len(page2["results"]) == 1
    assert page2["results"][0]["annotation_id"] == "ann-mouse"
    assert page2["previous"]


def test_annotation_genes_browse(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes")
    assert r.status_code == 200
    body = r.json()
    assert body["annotation_id"] == "ann-human"
    assert [x["local_id"] for x in body["results"]] == [10, 20]
    for row in body["results"]:
        assert "xrefs" not in row
        assert "matched" not in row
        assert row["annotation_id"] == "ann-human"


def test_annotation_namespaces(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/namespaces")
    assert r.status_code == 200
    body = r.json()
    assert body["annotation_id"] == "ann-human"
    assert body["next"] is None
    assert body["previous"] is None
    by_ns = {row["namespace"]: row["accession_count"] for row in body["results"]}
    assert by_ns == {
        "alias": 1,
        "ensembl_transcript": 1,
        "go": 1,
        "symbol": 2,
    }
    assert body["total"] == 4
    assert [row["namespace"] for row in body["results"]] == sorted(by_ns)


def test_annotation_namespaces_missing_shard(client: TestClient) -> None:
    r = client.get("/annotations/missing/namespaces")
    assert r.status_code == 404
    assert r.json()["detail"]["code"] == "shard_not_found"


def test_annotation_genes_q_is_name_prefix_not_curie(client: TestClient) -> None:
    # Colon is part of the name prefix — no longer an xref CURIE filter.
    r = client.get(
        "/annotations/ann-human/genes",
        params={"q": "symbol:tp53"},
    )
    assert r.status_code == 200
    assert r.json()["results"] == []

    bad_prefix = client.get(
        "/annotations/ann-human/genes",
        params={"q": "notans:foo"},
    )
    assert bad_prefix.status_code == 200
    assert bad_prefix.json()["results"] == []


def test_annotation_genes_curies(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53"},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 1
    assert body["results"][0]["local_id"] == 10
    assert "matched" not in body["results"][0]
    assert "xrefs" not in body["results"][0]
    assert "errors" not in body


def test_annotation_genes_curies_match_any(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,symbol:brca1", "match": "any"},
    )
    assert r.status_code == 200
    assert [x["local_id"] for x in r.json()["results"]] == [10, 20]


def test_annotation_genes_curies_match_all(client: TestClient) -> None:
    both = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,alias:tp53", "match": "all"},
    )
    assert both.status_code == 200
    assert [x["local_id"] for x in both.json()["results"]] == [10]

    missing = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,symbol:brca1", "match": "all"},
    )
    assert missing.status_code == 200
    assert missing.json()["results"] == []


def test_annotation_genes_curies_unknown_prefix(client: TestClient) -> None:
    any_r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,notans:foo", "match": "any"},
    )
    assert any_r.status_code == 200
    body = any_r.json()
    assert [x["local_id"] for x in body["results"]] == [10]
    assert body["errors"] == [
        {
            "curie": "notans:foo",
            "code": "unknown_prefix",
            "message": "unknown CURIE prefix: notans",
        }
    ]

    all_r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53,notans:foo", "match": "all"},
    )
    assert all_r.status_code == 200
    all_body = all_r.json()
    assert all_body["results"] == []
    assert all_body["errors"][0]["code"] == "unknown_prefix"


def test_annotation_genes_curies_required(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": " , "},
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "curies_required"


def test_annotation_genes_too_many_curies(client: TestClient) -> None:
    curies = ",".join(f"symbol:g{i}" for i in range(21))
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": curies},
    )
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert detail["code"] == "too_many_curies"
    assert detail["max"] == 20


def test_annotation_genes_invalid_match(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"curies": "symbol:tp53", "match": "nope"},
    )
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "invalid_match"


def test_annotation_genes_q_and_curies(client: TestClient) -> None:
    hit = client.get(
        "/annotations/ann-human/genes",
        params={"q": "TP", "curies": "symbol:tp53"},
    )
    assert hit.status_code == 200
    assert [x["local_id"] for x in hit.json()["results"]] == [10]

    miss = client.get(
        "/annotations/ann-human/genes",
        params={"q": "BRCA", "curies": "symbol:tp53"},
    )
    assert miss.status_code == 200
    assert miss.json()["results"] == []


def test_annotation_genes_q_primary_name(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"q": "TP53"},
    )
    assert r.status_code == 200
    assert [x["local_id"] for x in r.json()["results"]] == [10]

    prefix = client.get(
        "/annotations/ann-human/genes",
        params={"q": "tp"},
    )
    assert prefix.status_code == 200
    assert [x["local_id"] for x in prefix.json()["results"]] == [10]


def test_annotation_genes_q_too_long(client: TestClient) -> None:
    r = client.get(
        "/annotations/ann-human/genes",
        params={"q": "a" * 65},
    )
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert detail["code"] == "q_too_long"
    assert detail["max"] == 64
    assert detail["received"] == 65


def test_annotation_genes_missing_shard(client: TestClient) -> None:
    r = client.get("/annotations/missing/genes")
    assert r.status_code == 404
    assert r.json()["detail"]["code"] == "shard_not_found"


def test_annotation_genes_pagination(client: TestClient) -> None:
    first = client.get(
        "/annotations/ann-human/genes",
        params={"limit": 1},
    ).json()
    assert len(first["results"]) == 1
    assert first["results"][0]["local_id"] == 10
    assert first["next"]
    assert first["previous"] is None

    second = client.get(
        "/annotations/ann-human/genes",
        params={"limit": 1, "next": first["next"]},
    ).json()
    assert len(second["results"]) == 1
    assert second["results"][0]["local_id"] == 20
    assert second["previous"]
    assert second["next"] is None

    back = client.get(
        "/annotations/ann-human/genes",
        params={"limit": 1, "previous": second["previous"]},
    ).json()
    assert back["results"][0]["local_id"] == 10


def test_gene_xrefs_list(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/10/xrefs")
    assert r.status_code == 200
    body = r.json()
    assert body["annotation_id"] == "ann-human"
    assert body["local_id"] == 10
    ns = {x["namespace"] for x in body["results"]}
    assert ns == {"alias", "ensembl_transcript", "go", "symbol"}
    # Ordered by namespace, accession
    assert [x["namespace"] for x in body["results"]] == [
        "alias",
        "ensembl_transcript",
        "go",
        "symbol",
    ]


def test_gene_xrefs_missing(client: TestClient) -> None:
    r = client.get("/annotations/ann-human/genes/999/xrefs")
    assert r.status_code == 404
    assert r.json()["detail"]["code"] == "gene_not_found"
    r2 = client.get("/annotations/missing/genes/10/xrefs")
    assert r2.status_code == 404
    assert r2.json()["detail"]["code"] == "shard_not_found"


def test_gene_xrefs_pagination(client: TestClient) -> None:
    first = client.get(
        "/annotations/ann-human/genes/10/xrefs",
        params={"limit": 1},
    ).json()
    assert len(first["results"]) == 1
    assert first["results"][0]["namespace"] == "alias"
    assert first["next"]
    assert first["previous"] is None

    second = client.get(
        "/annotations/ann-human/genes/10/xrefs",
        params={"limit": 1, "next": first["next"]},
    ).json()
    assert len(second["results"]) == 1
    assert second["results"][0]["namespace"] == "ensembl_transcript"
    assert second["previous"]

    back = client.get(
        "/annotations/ann-human/genes/10/xrefs",
        params={"limit": 1, "previous": second["previous"]},
    ).json()
    assert back["results"][0]["namespace"] == "alias"
