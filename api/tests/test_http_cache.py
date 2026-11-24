"""ETag/Cache-Control by data version, and the gffread child limit (API-LOCKED-DECISIONS.md, D9)."""
from __future__ import annotations

import subprocess
import sys

import helpers.gffread_runner as gr
import helpers.presence as presence_mod
import settings as settings_mod
from tests.test_presence_api import presence_client  # noqa: F401


def test_etag_cache_control_and_304(presence_client):  # noqa: F811
    client, _ = presence_client
    r = client.get("/hits", params={"names": "tp53"})
    assert r.status_code == 200
    assert r.headers["etag"] == 'W/"abc123"'
    assert r.headers["cache-control"] == "public, max-age=300"
    again = client.get("/hits", params={"names": "tp53"}, headers={"If-None-Match": r.headers["etag"]})
    assert again.status_code == 304 and again.content == b""
    assert again.headers["etag"] == r.headers["etag"]
    other = client.get("/hits", params={"names": "tp53"}, headers={"If-None-Match": 'W/"old"'})
    assert other.status_code == 200


def test_no_cache_headers_on_errors_health_and_other_methods(presence_client):  # noqa: F811
    client, _ = presence_client
    bad = client.get("/hits", params={"names": "tp53", "shape": "nope"})
    assert bad.status_code >= 400 and "etag" not in bad.headers
    assert "etag" not in client.get("/health").headers
    assert client.post("/hits", json={"names": ["tp53"]}).headers.get("etag") is None


def test_version_change_during_request_is_not_cached(presence_client, monkeypatch):  # noqa: F811
    client, _ = presence_client
    client.get("/health")
    seen = iter(["abc123", "newer"])
    import helpers.http_cache as hc
    monkeypatch.setattr(hc, "current_data_version", lambda: next(seen, "newer"))
    r = client.get("/hits", params={"names": "tp53"})
    assert r.headers["cache-control"] == "no-store" and "etag" not in r.headers


def test_max_age_zero_disables_cache_control(presence_client, monkeypatch):  # noqa: F811
    client, _ = presence_client
    monkeypatch.setattr(settings_mod.settings, "HTTP_CACHE_MAX_AGE_S", 0)
    r = client.get("/hits", params={"names": "tp53"})
    assert "cache-control" not in r.headers and r.headers["etag"]


def test_gffread_child_is_limited(monkeypatch):
    monkeypatch.setattr(settings_mod.settings, "GFFREAD_MAX_DATA_MB", 64)
    ok = subprocess.run([sys.executable, "-c", "x = bytearray(8 * 1024 * 1024)"],
                        preexec_fn=gr.limit_child, capture_output=True)
    assert ok.returncode == 0
    big = subprocess.run([sys.executable, "-c", "x = bytearray(400 * 1024 * 1024)"],
                         preexec_fn=gr.limit_child, capture_output=True)
    assert big.returncode != 0 and b"MemoryError" in big.stderr
