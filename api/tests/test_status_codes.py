"""One status for every rejected input (422) and the shape/format rules of /hits."""
from __future__ import annotations

import re
from pathlib import Path

from tests.test_presence_api import presence_client  # noqa: F401


def test_no_route_file_answers_400():
    for path in (Path(__file__).parent.parent / "routers").glob("*.py"):
        assert "status_code=400" not in path.read_text(), path.name
        assert not re.search(r"\bstatus=400\b", path.read_text()), path.name


def test_rejected_input_is_422_everywhere(presence_client):  # noqa: F811
    client, _ = presence_client
    cases = [
        ("/hits", {"names": "tp53", "alias": "nope"}),
        ("/hits", {"names": "tp53", "shape": "nope"}),
        ("/hits", {"names": "tp53", "format": "xml"}),
        ("/hits", {"names": "tp53", "min_matched": 9}),
        ("/hits", {"names": "se*"}),                       # prefix too short
        ("/names", {"prefix": "se"}),
        ("/names", {"prefix": "tp53", "sources": "Nowhere"}),
        ("/annotations/ann-human/genes", {"offset": 99_999_999}),
    ]
    for path, params in cases:
        r = client.get(path, params=params)
        assert r.status_code == 422, (path, params, r.status_code, r.text[:200])
        detail = r.json()["detail"]
        assert detail["code"] and detail["message"], (path, params)


def test_tsv_with_compact_shape_is_refused_not_ignored(presence_client):  # noqa: F811
    client, _ = presence_client
    r = client.get("/hits", params={"names": "tp53", "format": "tsv", "shape": "compact"})
    assert r.status_code == 422 and r.json()["detail"]["code"] == "invalid_combination"
    ok = client.get("/hits", params={"names": "tp53", "format": "tsv", "shape": "records"})
    assert ok.status_code == 200 and ok.text.startswith("#format=tsv")
    plain = client.get("/hits", params={"names": "tp53", "format": "tsv"})
    assert plain.status_code == 200


def test_shape_has_no_default_and_is_echoed(presence_client):  # noqa: F811
    client, _ = presence_client
    body = client.get("/hits", params={"names": "tp53"}).json()
    assert body["shape"] == "compact" and "rows" in body
    rec = client.get("/hits", params={"names": "tp53", "shape": "records"}).json()
    assert rec["shape"] == "records" and "results" in rec
    spec = client.get("/openapi.json").json()
    shape = next(p for p in spec["paths"]["/hits"]["get"]["parameters"] if p["name"] == "shape")
    assert "default" not in shape["schema"] and shape["schema"].get("anyOf")


def test_a_list_answer_has_rows_or_results_never_both(presence_client):  # noqa: F811
    client, _ = presence_client
    compact = client.get("/hits", params={"names": "tp53"}).json()
    assert "rows" in compact and "results" not in compact and "fields" in compact
    records = client.get("/hits", params={"names": "tp53", "shape": "records"}).json()
    assert "results" in records and "rows" not in records
    names = client.get("/names", params={"prefix": "tp5"}).json()
    assert "results" in names and "rows" not in names
