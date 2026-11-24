"""Concurrency gate and time budget of /annotations/*."""
from __future__ import annotations

import sqlite3
import time

import pytest
from fastapi.testclient import TestClient

import helpers.concurrency as conc
import services.annotations as ann_mod
import settings as settings_mod
from main import create_app


def test_gate_refuses_when_all_slots_taken(monkeypatch):
    sem = conc.threading.Semaphore(1)
    monkeypatch.setattr(conc, "_ANNOTATIONS_SEMAPHORE", sem)
    with TestClient(create_app()) as client:
        sem.acquire()
        r = client.get("/annotations/abc/namespaces")
        sem.release()
        assert r.status_code == 503
        assert r.headers["retry-after"] == "1"
        assert r.json()["detail"]["code"] == "annotations_busy"
        # slot is free again: the route is reached (shard is simply missing here -> 404/503 not busy)
        r2 = client.get("/annotations/abc/namespaces")
        assert r2.json()["detail"]["code"] != "annotations_busy"


def test_slot_is_released_after_request(monkeypatch):
    sem = conc.threading.Semaphore(1)
    monkeypatch.setattr(conc, "_ANNOTATIONS_SEMAPHORE", sem)
    with TestClient(create_app()) as client:
        for _ in range(3):
            client.get("/annotations/abc/namespaces")
    assert sem.acquire(blocking=False)


def test_time_budget_interrupts_a_long_query(tmp_path, monkeypatch):
    db = tmp_path / "genes.sqlite"
    c = sqlite3.connect(db)
    c.execute("CREATE TABLE t(x)")
    c.commit()
    c.close()
    monkeypatch.setattr(ann_mod, "per_gff_root", lambda base=None: tmp_path)
    monkeypatch.setattr(ann_mod, "shard_genes_path", lambda *_a, **_k: db)
    monkeypatch.setattr(settings_mod.settings, "ANNOTATIONS_QUERY_BUDGET_S", 0.2)
    conn = ann_mod.open_shard("x")
    t0 = time.monotonic()
    with pytest.raises(sqlite3.OperationalError, match="interrupted"):
        conn.execute("WITH RECURSIVE r(n) AS (SELECT 1 UNION ALL SELECT n+1 FROM r) "
                     "SELECT count(*) FROM r").fetchone()
    assert time.monotonic() - t0 < 5
    conn.close()


def test_threadpool_is_capped():
    import anyio.to_thread
    with TestClient(create_app()) as client:
        client.get("/health")
        # the limiter belongs to the app's event loop: read it from inside it
        async def tokens():
            return anyio.to_thread.current_default_thread_limiter().total_tokens
        assert client.portal.call(tokens) == settings_mod.settings.THREADPOOL_SIZE == 12
