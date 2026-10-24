"""Tests for GET /names prefix discovery."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from helpers.pagination import MAX_PREFIX_EXPAND
from services.names import list_names_by_prefix


def test_names_prefix(client: TestClient) -> None:
    r = client.get("/names", params={"prefix": "gpx4"})
    assert r.status_code == 200
    body = r.json()
    assert body["prefix"] == "gpx4"
    assert body["include_taxids"] == []
    assert body["exclude_taxids"] == []
    assert body["total"] == 3
    names = [row["name"] for row in body["results"]]
    assert names == ["gpx4", "gpx4a", "gpx4b"]
    by_name = {row["name"]: row for row in body["results"]}
    assert by_name["gpx4"]["n_annotations"] == 2
    assert by_name["gpx4"]["category"] == "exact"
    assert by_name["gpx4a"]["category"] == "exact"
    assert by_name["gpx4b"]["n_annotations"] == 1


def test_names_categories(client: TestClient) -> None:
    # tp53 is exact; use a synthetic check via classify on known shapes
    from helpers.name_category import classify_name

    assert classify_name("selenop") == "exact"
    assert classify_name("gpx4-1") == "numbered"
    assert classify_name("gpx4.a") == "homeolog"
    assert classify_name("gpx4.ps") == "pseudogene"
    assert classify_name("gpx4_ps") == "pseudogene"
    assert classify_name("predicted_gpx4") == "predicted"
    assert classify_name("gpx4_like") == "other"


def test_names_include_taxids(client: TestClient) -> None:
    r = client.get(
        "/names",
        params={"prefix": "gpx4", "include_taxids": "9606"},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["include_taxids"] == [9606]
    # human has gpx4 + gpx4a; mouse has gpx4 + gpx4b
    names = {row["name"] for row in body["results"]}
    assert names == {"gpx4", "gpx4a"}
    by_name = {row["name"]: row for row in body["results"]}
    assert by_name["gpx4"]["n_annotations"] == 1


def test_names_exclude_taxids(client: TestClient) -> None:
    r = client.get(
        "/names",
        params={"prefix": "gpx4", "exclude_taxids": "9606"},
    )
    assert r.status_code == 200
    body = r.json()
    names = {row["name"] for row in body["results"]}
    assert names == {"gpx4", "gpx4b"}


def test_names_validation(client: TestClient) -> None:
    assert client.get("/names").status_code == 422

    short = client.get("/names", params={"prefix": "ab"})
    assert short.status_code == 422
    assert short.json()["detail"]["code"] == "prefix_too_short"

    star = client.get("/names", params={"prefix": "gpx4*"})
    assert star.status_code == 422
    assert star.json()["detail"]["code"] == "invalid_prefix"

    bad_tax = client.get(
        "/names", params={"prefix": "gpx4", "include_taxids": "nope"}
    )
    assert bad_tax.status_code == 422
    assert bad_tax.json()["detail"]["code"] == "invalid_taxid"


def test_names_prefix_too_broad_http(
    client: TestClient, monkeypatch: pytest.MonkeyPatch
) -> None:
    """HTTP path returns 422 when expand cap is exceeded."""
    import helpers.names as names_helpers
    import services.names as names_svc

    monkeypatch.setattr(names_helpers, "MAX_PREFIX_EXPAND", 2)
    monkeypatch.setattr(names_svc, "MAX_PREFIX_EXPAND", 2)
    r = client.get("/names", params={"prefix": "gpx4"})
    assert r.status_code == 422
    assert r.json()["detail"]["code"] == "prefix_too_broad"


def test_names_prefix_too_broad(tmp_path: Path) -> None:
    db = tmp_path / "broad.sqlite"
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    conn.executescript(
        """
        CREATE TABLE xref_file (
            namespace TEXT, accession TEXT, annotation_key INT,
            PRIMARY KEY (namespace, accession, annotation_key)
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
            "INSERT INTO xref_file VALUES ('symbol', ?, 1)", (acc,)
        )
    conn.execute("INSERT INTO annotation_lineage VALUES (1, 9606)")
    conn.commit()
    with pytest.raises(Exception) as excinfo:
        list_names_by_prefix(conn, "zzz")
    conn.close()
    from helpers.cursor import QueryError

    assert isinstance(excinfo.value, QueryError)
    assert excinfo.value.code == "prefix_too_broad"
