"""Unexpected failures answer with the documented error shape and leave a log line."""
from __future__ import annotations

import logging
import sqlite3

import pytest
from fastapi.testclient import TestClient

from main import create_app


def _client(exc: Exception) -> TestClient:
    app = create_app()

    @app.get("/boom")
    def boom():
        raise exc

    return TestClient(app, raise_server_exceptions=False)


def test_an_unexpected_exception_is_a_coded_500(caplog):
    with caplog.at_level(logging.ERROR):
        r = _client(ValueError("secret path /data/x")).get("/boom")
    assert r.status_code == 500
    assert r.json() == {"detail": {"code": "internal_error", "message": "Internal Server Error"}}
    assert "/data/x" not in r.text                      # the cause stays in the log


def test_a_sqlite_failure_is_logged_and_coded(caplog):
    with caplog.at_level(logging.ERROR):
        r = _client(sqlite3.OperationalError("no such table: gene")).get("/boom")
    assert r.status_code == 500 and r.json()["detail"]["code"] == "internal_error"
    assert any("no such table" in rec.getMessage() for rec in caplog.records)


@pytest.mark.parametrize("message", ["disk I/O error", "unable to open database file"])
def test_a_storage_hiccup_is_a_retryable_503(message):
    r = _client(sqlite3.OperationalError(message)).get("/boom")
    assert r.status_code == 503 and r.headers["retry-after"] == "5"
    assert r.json()["detail"]["code"] == "index_unavailable"


def test_the_query_budget_is_still_a_503_query_too_slow():
    r = _client(sqlite3.OperationalError("interrupted")).get("/boom")
    assert r.status_code == 503 and r.json()["detail"]["code"] == "query_too_slow"
