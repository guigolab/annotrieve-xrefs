"""API integration tests for the reshaped v1 surface."""
from __future__ import annotations

import sqlite3

import pytest
from fastapi.testclient import TestClient


def test_health(client: TestClient) -> None:
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_namespaces_allowlist_only(client: TestClient) -> None:
    r = client.get("/namespaces")
    assert r.status_code == 200
    body = r.json()
    assert "offset" in body
    assert body["total"] == len(body["results"])
    ns = {row["namespace"] for row in body["results"]}
    assert ns == {"alias", "go", "interpro", "symbol"}
    for row in body["results"]:
        assert "labels" in row
        assert "match_count" not in row
    go = next(row for row in body["results"] if row["namespace"] == "go")
    assert go["labels"] is True
    assert go["name"] == "Gene Ontology"


def test_namespaces_no_has_labels_param(client: TestClient) -> None:
    r = client.get("/namespaces", params={"has_labels": "true"})
    # Unused query params are ignored by FastAPI when not declared.
    assert r.status_code == 200
    assert all("match_count" not in row for row in r.json()["results"])


def test_namespace_item(client: TestClient) -> None:
    r = client.get("/namespaces/symbol")
    assert r.status_code == 200
    assert r.json() == {
        "namespace": "symbol",
        "name": "Gene symbol",
        "accession_count": 1,
        "labels": False,
    }
    go = client.get("/namespaces/go")
    assert go.json()["labels"] is True
    missing = client.get("/namespaces/nope")
    assert missing.status_code == 404
    demoted = client.get("/namespaces/ensembl")
    assert demoted.status_code == 404


def test_matches_requires_q(client: TestClient) -> None:
    r = client.get("/matches")
    assert r.status_code == 422


def test_matches_prose(client: TestClient) -> None:
    r = client.get("/matches", params={"q": "biological"})
    assert r.status_code == 200
    body = r.json()
    assert body["total"] >= 1
    row = body["results"][0]
    assert row["namespace"] == "go"
    assert row["match_count"] >= 1
    assert "labels" in row


def test_matches_symbol_token(client: TestClient) -> None:
    r = client.get("/matches", params={"q": "tp53"})
    assert r.status_code == 200
    by_ns = {row["namespace"]: row["match_count"] for row in r.json()["results"]}
    assert by_ns.get("symbol", 0) >= 1
    assert by_ns.get("alias", 0) >= 1


def test_matches_short_q(client: TestClient) -> None:
    r = client.get("/matches", params={"q": "tp"})
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "q_too_short"


def test_matches_curie(client: TestClient) -> None:
    r = client.get("/matches", params={"q": "symbol:tp53"})
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["results"][0]["namespace"] == "symbol"
    assert body["results"][0]["match_count"] == 1


def test_terms_go_browse(client: TestClient) -> None:
    r = client.get("/namespaces/go/terms")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["offset"] == 0
    assert body["results"][0]["accession"] == "GO:0008150"
    assert body["results"][0]["label"] == "biological_process"
    assert body["results"][0]["curie"] == "go:GO:0008150"


def test_terms_go_search(client: TestClient) -> None:
    r = client.get("/namespaces/go/terms", params={"q": "biological"})
    assert r.status_code == 200
    assert r.json()["total"] == 1
    assert r.json()["results"][0]["n_genes"] == 2


def test_terms_symbol_requires_q(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/terms")
    assert r.status_code == 400
    assert r.json()["detail"]["code"] == "q_required"


def test_terms_symbol_prefix(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/terms", params={"q": "tp53"})
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["results"][0]["accession"] == "tp53"


def test_term_detail(client: TestClient) -> None:
    r = client.get("/namespaces/symbol/terms/tp53")
    assert r.status_code == 200
    assert r.json()["curie"] == "symbol:tp53"
    assert r.json()["n_annotations"] == 2
    folded = client.get("/namespaces/symbol/terms/TP53")
    assert folded.status_code == 200
    go = client.get("/namespaces/go/terms/8150")
    assert go.status_code == 200
    assert go.json()["accession"] == "GO:0008150"
    missing = client.get("/namespaces/symbol/terms/nope")
    assert missing.status_code == 404
    demoted = client.get("/namespaces/ensembl/terms/ENSG1")
    assert demoted.status_code == 404


def test_old_labels_and_accessions_gone(client: TestClient) -> None:
    assert client.get("/namespaces/go/labels").status_code == 404
    assert client.get("/namespaces/symbol/accessions").status_code == 404
    assert client.get("/hits").status_code == 404
    assert client.get("/annotations/index").status_code == 404


def test_annotations_catalogue(client: TestClient) -> None:
    r = client.get("/annotations")
    assert r.status_code == 200
    body = r.json()
    assert body["n_without_tier_a"] == 0
    assert body["fields"][0] == "annotation_id"
    assert len(body["rows"]) == 2
    assert body["rows"][0][1] == 1


def test_annotations_presence_search(client: TestClient) -> None:
    r = client.get("/annotations", params={"curies": "symbol:tp53"})
    assert r.status_code == 200
    body = r.json()
    assert body["fields"] == [
        "annotation_id",
        "annotation_key",
        "taxid",
        "source",
    ]
    keys = {row[1] for row in body["rows"]}
    assert keys == {1, 2}


def test_annotations_presence_post(client: TestClient) -> None:
    r = client.post(
        "/annotations",
        json={"curies": ["symbol:tp53"], "taxid": 9606},
    )
    assert r.status_code == 200
    assert len(r.json()["rows"]) == 1
    assert r.json()["taxid"] == 9606


def test_annotations_presence_match_all(client: TestClient) -> None:
    r = client.get(
        "/annotations",
        params={"curies": "symbol:tp53,alias:tp53", "match": "all"},
    )
    assert r.status_code == 200
    # Only human has both symbol and alias tp53 in xref_file.
    assert [row[0] for row in r.json()["rows"]] == ["ann-human"]


def test_annotations_unknown_prefix(client: TestClient) -> None:
    r = client.get(
        "/annotations",
        params={"curies": "ensembl:ENSG1", "match": "any"},
    )
    assert r.status_code == 200
    assert r.json()["rows"] == []
    assert r.json()["errors"][0]["code"] == "unknown_prefix"


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


def test_annotations_busy_returns_503(
    client: TestClient, monkeypatch: pytest.MonkeyPatch
) -> None:
    from helpers.concurrency import BusyError

    class _Busy:
        def __enter__(self):
            raise BusyError()

        def __exit__(self, *args):
            return False

    monkeypatch.setattr("routers.annotations.hits_slot", lambda: _Busy())
    r = client.get("/annotations", params={"curies": "symbol:tp53"})
    assert r.status_code == 503
    assert r.headers.get("Retry-After") == "1"


def test_matches_and_terms_accession_id(client: TestClient) -> None:
    """Labeled ns: paste accession id when FTS text does not contain it."""
    matches = client.get("/matches", params={"q": "IPR002117"})
    assert matches.status_code == 200
    by_ns = {
        row["namespace"]: row["match_count"]
        for row in matches.json()["results"]
    }
    assert by_ns.get("interpro", 0) >= 1

    terms = client.get(
        "/namespaces/interpro/terms", params={"q": "IPR002117"}
    )
    assert terms.status_code == 200
    body = terms.json()
    assert body["total"] >= 1
    assert body["results"][0]["accession"] == "IPR002117"


def _seed_pfam_versioned(corpus_dir) -> None:
    db = corpus_dir / "gene_corpus.sqlite"
    conn = sqlite3.connect(str(db))
    conn.execute("INSERT INTO namespace_stats VALUES ('pfam', 1)")
    conn.execute(
        "INSERT INTO xref_meta VALUES ('pfam', 'PF00001.15', 3, 10)"
    )
    conn.execute(
        "INSERT INTO xref_label VALUES ('pfam', 'PF00001', '7tm_1')"
    )
    try:
        conn.execute(
            "INSERT INTO xref_label_fts (namespace, accession, text) "
            "VALUES ('pfam', 'PF00001', '7tm_1')"
        )
    except sqlite3.OperationalError:
        pass
    conn.commit()
    conn.close()


def test_terms_pfam_versioned_meta_counts(
    client: TestClient, corpus_dir
) -> None:
    _seed_pfam_versioned(corpus_dir)

    r = client.get("/namespaces/pfam/terms")
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["results"][0]["accession"] == "PF00001"
    assert body["results"][0]["n_annotations"] == 3
    assert body["results"][0]["n_genes"] == 10


def test_pfam_free_text_versioned_and_unversioned(
    client: TestClient, corpus_dir
) -> None:
    _seed_pfam_versioned(corpus_dir)

    for q in ("PF00001", "PF00001.15"):
        matches = client.get("/matches", params={"q": q})
        assert matches.status_code == 200
        by_ns = {
            row["namespace"]: row["match_count"]
            for row in matches.json()["results"]
        }
        assert by_ns.get("pfam", 0) >= 1, q

        terms = client.get("/namespaces/pfam/terms", params={"q": q})
        assert terms.status_code == 200
        body = terms.json()
        assert body["total"] >= 1, q
        assert body["results"][0]["accession"] == "PF00001"


def test_pfam_exact_term_requires_versioned_meta(
    client: TestClient, corpus_dir
) -> None:
    _seed_pfam_versioned(corpus_dir)

    missing = client.get("/namespaces/pfam/terms/PF00001")
    assert missing.status_code == 404

    found = client.get("/namespaces/pfam/terms/PF00001.15")
    assert found.status_code == 200
    assert found.json()["accession"] == "PF00001.15"
    assert found.json()["n_annotations"] == 3
