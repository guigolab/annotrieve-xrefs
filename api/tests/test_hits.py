"""Tests for GET|POST /hits gene-name presence search."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from services.hits import list_name_hits


def test_hits_exact_name(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "TP53"})
    assert r.status_code == 200
    body = r.json()
    assert body["fields"] == [
        "annotation_id",
        "taxid",
        "source",
        "matched",
        "release_date",
        "assembly_name",
    ]
    assert "keys" not in body
    assert body["sources"]["3"] == {"label": "RefSeq"}
    assert body["include_taxids"] == []
    assert body["exclude_taxids"] == []
    assert len(body["terms"]) == 1
    term = body["terms"][0]
    assert term["input"] == "TP53"
    assert term["name"] == "tp53"
    assert set(term["namespaces"]) == {"alias", "symbol"}
    assert term["n_hits"] == 2
    assert "expanded" not in term

    assert body["total"] == 2
    assert len(body["rows"]) == 2
    human = next(row for row in body["rows"] if row[0] == "ann-human")
    assert human[1] == 9606
    assert human[2] == 3  # RefSeq
    assert human[3] == [0]  # term indices
    assert human[4] == "2024-01-02T00:00:00"
    assert human[5] == "GRCh38.p14"
    mouse = next(row for row in body["rows"] if row[0] == "ann-mouse")
    assert mouse[3] == [0]
    assert mouse[4] is None
    assert mouse[5] is None


def test_hits_rejects_wildcard(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "gpx4*"})
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "wildcard_not_supported"


def test_hits_multi_term_union(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "selenop,sepp1"})
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 1
    row = body["rows"][0]
    assert row[0] == "ann-human"
    assert set(row[3]) == {0, 1}


def test_hits_unknown_name(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "tp53,nope"})
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 2
    assert body["terms"][1]["name"] == "nope"
    assert body["terms"][1]["n_hits"] == 0
    assert body["terms"][1]["namespaces"] == []


def test_hits_dedupe_casefold(client: TestClient) -> None:
    r = client.get("/hits", params={"names": "TP53,tp53"})
    assert r.status_code == 200
    body = r.json()
    assert len(body["terms"]) == 1
    assert body["terms"][0]["input"] == "TP53"


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
    assert set(body["terms"][0]["namespaces"]) == {"alias", "symbol"}

    no_human = client.get(
        "/hits",
        params={"names": "tp53", "exclude_taxids": "9606"},
    )
    assert no_human.status_code == 200
    body = no_human.json()
    assert body["exclude_taxids"] == [9606]
    assert body["total"] == 1
    assert body["rows"][0][0] == "ann-mouse"
    assert body["terms"][0]["namespaces"] == ["symbol"]


def test_hits_post_and_multi_term_matched(client: TestClient) -> None:
    r = client.post(
        "/hits",
        json={"names": ["SELENOP", "SEPP1", "gpx4"]},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["terms"]) == 3
    assert "expanded" not in body["terms"][2]
    human = next(row for row in body["rows"] if row[0] == "ann-human")
    assert set(human[3]) == {0, 1, 2}


def test_hits_validation_errors(client: TestClient) -> None:
    assert client.get("/hits").status_code == 422  # missing names

    empty = client.get("/hits", params={"names": " , "})
    assert empty.status_code == 422
    assert empty.json()["detail"]["code"] == "names_required"

    bad_tax = client.get(
        "/hits", params={"names": "tp53", "include_taxids": "nope"}
    )
    assert bad_tax.status_code == 422
    assert bad_tax.json()["detail"]["code"] == "invalid_taxid"


def test_hits_too_many_names(client: TestClient) -> None:
    from helpers.pagination import MAX_NAMES_GET

    names = ",".join(f"gene{i}" for i in range(MAX_NAMES_GET + 1))
    r = client.get("/hits", params={"names": names})
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "too_many_names"
    assert "n_matched" not in (
        client.get("/hits", params={"names": "tp53"}).json()["fields"]
    )


def test_hits_post_too_many_names(client: TestClient) -> None:
    """POST shares too_many_names with GET (no pydantic max_length shape)."""
    from helpers.pagination import MAX_NAMES_POST

    names = [f"gene{i}" for i in range(MAX_NAMES_POST + 1)]
    r = client.post("/hits", json={"names": names})
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "too_many_names"


def test_hits_post_casefold_dedupe_under_cap(client: TestClient) -> None:
    from helpers.pagination import MAX_NAMES_POST

    names: list[str] = []
    for i in range(MAX_NAMES_POST):
        names.append(f"gene{i}")
        names.append(f"GENE{i}")
    r = client.post("/hits", json={"names": names})
    assert r.status_code == 200
    assert len(r.json()["terms"]) == MAX_NAMES_POST


def test_hits_casefold_dedupe_under_cap(client: TestClient) -> None:
    """Casefolded duplicates do not count toward too_many_names."""
    from helpers.pagination import MAX_NAMES_GET

    parts: list[str] = []
    for i in range(MAX_NAMES_GET):
        parts.append(f"gene{i}")
        parts.append(f"GENE{i}")
    r = client.get("/hits", params={"names": ",".join(parts)})
    assert r.status_code == 200
    body = r.json()
    assert len(body["terms"]) == MAX_NAMES_GET
    assert "keys" not in body
    assert "n_matched" not in body["fields"]


def test_hits_orphan_xref_keys_ignored(tmp_path: Path) -> None:
    """xref_file keys missing from annotation do not inflate n_hits."""
    db = tmp_path / "orphan.sqlite"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(
        """
        CREATE TABLE xref_file (
            namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)
        );
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT, taxid INT, assembly_accession TEXT,
            organism_name TEXT, source_database TEXT,
            source_provider TEXT, profile_id TEXT, gff_path TEXT,
            release_date TEXT, assembly_name TEXT
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
        (1, 'ann-real', 9606, 'GCF_x', 'Homo', 'RefSeq', 'NCBI', 'refseq',
         '/t', NULL, NULL)
        """
    )
    conn.execute("INSERT INTO annotation_lineage VALUES (1, 9606)")
    conn.execute("INSERT INTO xref_file VALUES ('symbol', 'tp53', 1)")
    conn.execute("INSERT INTO xref_file VALUES ('symbol', 'tp53', 999)")  # orphan
    conn.commit()
    body = list_name_hits(conn, ["TP53"])
    conn.close()
    assert body["total"] == 1
    assert body["terms"][0]["n_hits"] == 1
    assert body["rows"][0][0] == "ann-real"


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


def test_hits_chunked_postings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Many accessions must not blow SQLite bind limits."""
    import services.hits as hits_mod

    monkeypatch.setattr(hits_mod, "_SQLITE_MAX_VARS", 6)

    db = tmp_path / "chunk.sqlite"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(
        """
        CREATE TABLE xref_file (
            namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)
        );
        CREATE TABLE annotation (
            annotation_key INTEGER PRIMARY KEY,
            annotation_id TEXT, taxid INT, assembly_accession TEXT,
            organism_name TEXT, source_database TEXT,
            source_provider TEXT, profile_id TEXT, gff_path TEXT,
            release_date TEXT, assembly_name TEXT
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
        (1, 'ann-a', 9606, 'GCF_x', 'Homo', 'RefSeq', 'NCBI', 'refseq',
         '/t', NULL, NULL)
        """
    )
    conn.execute("INSERT INTO annotation_lineage VALUES (1, 9606)")
    names = []
    for i in range(20):
        acc = f"gene{i:02d}"
        names.append(acc)
        conn.execute(
            "INSERT INTO xref_file VALUES ('symbol', ?, 1)", (acc,)
        )
    conn.commit()
    body = list_name_hits(conn, names, max_names=50)
    conn.close()
    assert body["total"] == 1
    assert body["terms"][0]["n_hits"] == 1
    assert len(body["terms"]) == 20
    assert set(body["rows"][0][3]) == set(range(20))
