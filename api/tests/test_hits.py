"""Tests for GET|POST /hits gene-name presence search."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from helpers.pagination import MAX_PREFIX_EXPAND
from services.hits import list_name_hits


def test_hits_exact_name(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "TP53"})
    assert r.status_code == 200
    body = r.json()
    assert body["fields"] == [
        "annotation_id",
        "taxid",
        "source",
        "n_matched",
        "matched",
    ]
    assert body["sources"]["3"] == {"label": "RefSeq"}
    assert body["include_taxids"] == []
    assert body["exclude_taxids"] == []
    assert len(body["terms"]) == 1
    term = body["terms"][0]
    assert term["input"] == "TP53"
    assert term["name"] == "tp53"
    assert set(term["matched_namespaces"]) == {"alias", "symbol"}
    assert term["n_hits"] == 2
    assert "expanded" not in term

    # keys: symbol+alias for tp53
    assert len(body["keys"]) == 2
    assert { (k["name"], k["ns"]) for k in body["keys"] } == {
        ("tp53", "symbol"),
        ("tp53", "alias"),
    }
    assert all(k["term"] == 0 for k in body["keys"])
    assert body["total"] == 2
    assert len(body["rows"]) == 2
    # Sorted by n_matched desc — human has both ns → more key hits
    human = next(row for row in body["rows"] if row[0] == "ann-human")
    assert human[1] == 9606
    assert human[2] == 3  # RefSeq
    assert human[3] == 1  # one term
    assert set(human[4]) == {0, 1}


def test_hits_prefix_expand(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "gpx4*"})
    assert r.status_code == 200
    body = r.json()
    term = body["terms"][0]
    assert term["name"] == "gpx4*"
    assert term["expanded"] == ["gpx4", "gpx4a", "gpx4b"]
    assert term["matched_namespaces"] == ["symbol"]
    assert term["n_hits"] == 2
    names = {k["name"] for k in body["keys"]}
    assert names == {"gpx4", "gpx4a", "gpx4b"}
    assert body["total"] == 2


def test_hits_match_all_and_any(client: TestClient) -> None:
    any_r = client.get(
        "/hits", params={"names": "selenop,sepp1", "match": "any"}
    )
    assert any_r.status_code == 200
    assert any_r.json()["total"] == 1

    all_r = client.get(
        "/hits", params={"names": "selenop,sepp1", "match": "all"}
    )
    assert all_r.status_code == 200
    body = all_r.json()
    assert body["total"] == 1
    row = body["rows"][0]
    assert row[0] == "ann-human"
    assert row[3] == 2  # both terms

    # missing name → match=all empty rows
    miss = client.get(
        "/hits", params={"names": "tp53,nope", "match": "all"}
    )
    assert miss.status_code == 200
    assert miss.json()["total"] == 0
    assert miss.json()["terms"][1]["n_hits"] == 0


def test_hits_include_exclude_taxids(client: TestClient) -> None:
    only_human = client.get(
        "/hits",
        params={"names": "tp53", "include_taxids": "9606"},
    )
    assert only_human.status_code == 200
    body = only_human.json()
    assert body["include_taxids"] == [9606]
    assert body["total"] == 1
    assert body["rows"][0][0] == "ann-human"
    assert body["terms"][0]["n_hits"] == 1

    no_human = client.get(
        "/hits",
        params={"names": "tp53", "exclude_taxids": "9606"},
    )
    assert no_human.status_code == 200
    body = no_human.json()
    assert body["exclude_taxids"] == [9606]
    assert body["total"] == 1
    assert body["rows"][0][0] == "ann-mouse"


def test_hits_post_and_multi_term_matched(client: TestClient) -> None:
    r = client.post(
        "/hits",
        json={
            "names": ["SELENOP", "SEPP1", "gpx4*"],
            "match": "any",
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["terms"]) == 3
    assert body["terms"][2]["expanded"] == ["gpx4", "gpx4a", "gpx4b"]
    # Human matches selenop + sepp1 + gpx4/gpx4a
    human = next(row for row in body["rows"] if row[0] == "ann-human")
    assert human[3] == 3
    key_terms = {body["keys"][i]["term"] for i in human[4]}
    assert key_terms == {0, 1, 2}


def test_hits_validation_errors(client: TestClient) -> None:
    assert client.get("/hits").status_code == 422  # missing names

    empty = client.get("/hits", params={"names": " , "})
    assert empty.status_code == 422
    assert empty.json()["detail"]["code"] == "names_required"

    short = client.get("/hits", params={"names": "ab*"})
    assert short.status_code == 422
    assert short.json()["detail"]["code"] == "prefix_too_short"

    bad_tax = client.get(
        "/hits", params={"names": "tp53", "include_taxids": "nope"}
    )
    assert bad_tax.status_code == 422
    assert bad_tax.json()["detail"]["code"] == "invalid_taxid"

    bad_match = client.get(
        "/hits", params={"names": "tp53", "match": "xor"}
    )
    assert bad_match.status_code == 422
    assert bad_match.json()["detail"]["code"] == "invalid_match"


def test_hits_tax_filter_strips_empty_namespaces(client: TestClient) -> None:
    """alias:tp53 is human-only; mouse include must drop alias from namespaces."""
    r = client.get(
        "/hits",
        params={"names": "tp53", "include_taxids": "10090"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    assert body["rows"][0][0] == "ann-mouse"
    term = body["terms"][0]
    assert term["matched_namespaces"] == ["symbol"]
    assert term["n_hits"] == 1
    # keys still list both ns (xref_file presence); only symbol contributes hits
    assert { (k["name"], k["ns"]) for k in body["keys"] } == {
        ("tp53", "symbol"),
        ("tp53", "alias"),
    }


def test_hits_too_many_taxids(client: TestClient) -> None:
    from helpers.pagination import MAX_TAXIDS_FILTER

    taxids = ",".join(str(i) for i in range(1, MAX_TAXIDS_FILTER + 2))
    r = client.get(
        "/hits",
        params={"names": "tp53", "include_taxids": taxids},
    )
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "too_many_taxids"

    r2 = client.post(
        "/hits",
        json={
            "names": ["tp53"],
            "exclude_taxids": list(range(1, MAX_TAXIDS_FILTER + 2)),
        },
    )
    assert r2.status_code == 422
    assert r2.json()["detail"]["code"] == "too_many_taxids"


def test_hits_post_int_taxids(client: TestClient) -> None:
    r = client.post(
        "/hits",
        json={
            "names": ["tp53"],
            "include_taxids": [9606],
            "exclude_taxids": [10090],
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert body["include_taxids"] == [9606]
    assert body["exclude_taxids"] == [10090]
    assert body["total"] == 1
    assert body["rows"][0][0] == "ann-human"


def test_hits_chunked_postings(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Many (ns, acc) pairs must not blow SQLite bind limits."""
    import services.hits as hits_mod

    monkeypatch.setattr(hits_mod, "_SQLITE_MAX_VARS", 6)

    db = tmp_path / "chunk.sqlite"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(
        """
        CREATE TABLE xref_meta (
            namespace TEXT, accession TEXT,
            n_annotations INT, n_loci INT,
            PRIMARY KEY (namespace, accession)
        );
        CREATE TABLE xref_file (
            namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)
        );
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT, taxid INT, assembly_accession TEXT,
            organism_name TEXT, source_database TEXT,
            source_provider TEXT, profile_id TEXT, gff_path TEXT
        );
        CREATE TABLE annotation_lineage (
            annotation_key INT, taxid INT,
            PRIMARY KEY (annotation_key, taxid)
        );
        """
    )
    conn.execute(
        """
        INSERT INTO annotation VALUES
        (1, 'ann-a', 9606, 'GCF_x', 'Homo', 'RefSeq', 'NCBI', 'refseq', '/t')
        """
    )
    conn.execute("INSERT INTO annotation_lineage VALUES (1, 9606)")
    names = []
    for i in range(20):
        acc = f"gene{i:02d}"
        names.append(acc)
        conn.execute(
            "INSERT INTO xref_meta VALUES ('symbol', ?, 1, 1)", (acc,)
        )
        conn.execute(
            "INSERT INTO xref_file VALUES ('symbol', ?, 1)", (acc,)
        )
    conn.commit()
    body = list_name_hits(conn, names, max_names=50)
    conn.close()
    assert body["total"] == 1
    assert body["terms"][0]["n_hits"] == 1
    assert len(body["keys"]) == 20


def test_hits_prefix_too_broad(tmp_path: Path) -> None:
    db = tmp_path / "broad.sqlite"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(
        """
        CREATE TABLE xref_meta (
            namespace TEXT, accession TEXT,
            n_annotations INT, n_loci INT,
            PRIMARY KEY (namespace, accession)
        );
        CREATE TABLE xref_file (
            namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)
        );
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT, taxid INT, assembly_accession TEXT,
            organism_name TEXT, source_database TEXT,
            source_provider TEXT, profile_id TEXT, gff_path TEXT
        );
        CREATE TABLE annotation_lineage (
            annotation_key INT, taxid INT,
            PRIMARY KEY (annotation_key, taxid)
        );
        """
    )
    for i in range(MAX_PREFIX_EXPAND + 1):
        acc = f"zzz{i:04d}"
        conn.execute(
            "INSERT INTO xref_meta VALUES ('symbol', ?, 1, 1)", (acc,)
        )
    conn.commit()
    with pytest.raises(Exception) as excinfo:
        list_name_hits(conn, "zzz*")
    conn.close()
    from helpers.cursor import QueryError

    assert isinstance(excinfo.value, QueryError)
    assert excinfo.value.code == "prefix_too_broad"
