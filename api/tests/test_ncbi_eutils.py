"""NCBI efetch client tests (mocked transport)."""
from __future__ import annotations

import time

import httpx
import pytest

from helpers.ncbi_eutils import (
    NcbiError,
    NcbiEutilsClient,
    reset_ncbi_client,
    close_ncbi_client,
)
from helpers.windows import GenomicWindow


def _fasta(seq: str, header: str = ">NC_1.1") -> str:
    return f"{header}\n{seq}\n"


@pytest.fixture(autouse=True)
def _reset_client():
    close_ncbi_client()
    yield
    close_ncbi_client()


def test_rate_without_key_is_slower(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("NCBI_API_KEY", "")
    monkeypatch.setenv("NCBI_WORKERS", "1")
    import settings as settings_mod

    settings_mod.settings.NCBI_API_KEY = ""
    settings_mod.settings.NCBI_WORKERS = 1

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=_fasta("ACGT"))

    transport = httpx.MockTransport(handler)
    client = NcbiEutilsClient(
        rate_per_sec=2.5,
        cache_mb=0,
        transport=transport,
    )
    assert abs(client.rate_interval - 0.4) < 1e-9
    assert not client.has_api_key
    client.close()


def test_rate_with_key(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("NCBI_API_KEY", "demo-key")
    import settings as settings_mod

    settings_mod.settings.NCBI_API_KEY = "demo-key"
    settings_mod.settings.NCBI_WORKERS = 1

    transport = httpx.MockTransport(
        lambda r: httpx.Response(200, text=_fasta("ACGT"))
    )
    client = NcbiEutilsClient(rate_per_sec=9.0, cache_mb=0, transport=transport)
    assert abs(client.rate_interval - (1 / 9.0)) < 1e-9
    assert client.has_api_key
    client.close()


def test_default_rate_not_divided_by_workers(monkeypatch: pytest.MonkeyPatch):
    """NCBI_WORKERS must not shrink the shared process rate bucket."""
    monkeypatch.setenv("NCBI_API_KEY", "demo-key")
    import settings as settings_mod

    settings_mod.settings.NCBI_API_KEY = "demo-key"
    settings_mod.settings.NCBI_WORKERS = 4

    transport = httpx.MockTransport(
        lambda r: httpx.Response(200, text=_fasta("ACGT"))
    )
    client = NcbiEutilsClient(cache_mb=0, transport=transport)
    # Full ~9/s budget, not 9/4.
    assert abs(client.rate_interval - (1 / 9.0)) < 1e-9
    client.close()


def test_fetch_success_and_cache(monkeypatch: pytest.MonkeyPatch):
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        return httpx.Response(200, text=_fasta("ACGTACGT"))

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=1,
        transport=httpx.MockTransport(handler),
    )
    reset_ncbi_client(client)
    window = GenomicWindow("NC_000001.11", 1, 8)
    a = client.fetch_window(window)
    b = client.fetch_window(window)
    assert a.sequence == "ACGTACGT"
    assert b.sequence == a.sequence
    assert calls["n"] == 1


def test_fetch_truncation_rejected():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=_fasta("ACGT"))

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )
    with pytest.raises(NcbiError) as exc:
        client.fetch_window(GenomicWindow("NC_000001.11", 1, 10))
    assert exc.value.code == "ncbi_sequence_truncated"
    client.close()


def test_fetch_truncation_allowed():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=_fasta("ACGT"))

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )
    slice_ = client.fetch_window(
        GenomicWindow("NC_000001.11", 1, 10),
        allow_truncation=True,
    )
    assert slice_.truncated
    assert slice_.sequence == "ACGT"
    client.close()


def test_non_fasta_200_is_error():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text="<html>error</html>")

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )
    with pytest.raises(NcbiError) as exc:
        client.fetch_window(GenomicWindow("NC_000001.11", 1, 4))
    assert exc.value.code == "ncbi_error"
    client.close()


def test_spaced_out_failure_body_is_sequence_not_found():
    body = "Error: F a i l e d  t o  u n d e r s t a n d  i d : scaffold_01"

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=body)

    client = NcbiEutilsClient(rate_per_sec=1000, cache_mb=0, transport=httpx.MockTransport(handler))
    with pytest.raises(NcbiError) as exc:
        client.fetch_window(GenomicWindow("scaffold_01", 1, 4))
    assert exc.value.code == "ncbi_sequence_not_found"
    assert "F a i l e d" not in str(exc.value)
    client.close()


def test_404_maps_to_not_found():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404, text="not found")

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )
    with pytest.raises(NcbiError) as exc:
        client.fetch_window(GenomicWindow("NC_000001.11", 1, 4))
    assert exc.value.code == "ncbi_sequence_not_found"
    client.close()


def test_429_retries_then_rate_limited():
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        return httpx.Response(429, headers={"Retry-After": "0"}, text="slow down")

    client = NcbiEutilsClient(
        rate_per_sec=1000,
        cache_mb=0,
        transport=httpx.MockTransport(handler),
    )
    with pytest.raises(NcbiError) as exc:
        client.fetch_window(GenomicWindow("NC_000001.11", 1, 4))
    assert exc.value.code == "ncbi_rate_limited"
    assert exc.value.status == 503
    assert calls["n"] == 3
    client.close()


# ---- one NCBI budget for all processes -------------------------------------------------------

def _worker(path: str, n: int, rate: float, q) -> None:
    from helpers.ncbi_eutils import SharedRateLimiter
    lim = SharedRateLimiter(rate, path)
    stamps = []
    for _ in range(n):
        lim.acquire()
        stamps.append(time.time())
    q.put(stamps)


def test_the_budget_is_shared_between_processes(tmp_path):
    """Four processes at 20 req/s each must together stay at 20 req/s, not 80."""
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    q = ctx.Queue()
    procs = [ctx.Process(target=_worker, args=(str(tmp_path / "rate"), 8, 20.0, q)) for _ in range(4)]
    for p in procs:
        p.start()
    stamps = sorted(t for _ in procs for t in q.get(timeout=30))
    for p in procs:
        p.join()
    assert len(stamps) == 32
    span = stamps[-1] - stamps[0]
    assert span >= 31 / 20.0 * 0.9, span                      # 31 gaps of 50 ms
    gaps = [b - a for a, b in zip(stamps, stamps[1:])]
    assert min(gaps) > 0.02                                    # never a burst


def test_threads_in_one_process_are_spaced_too(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from helpers.ncbi_eutils import SharedRateLimiter
    lim = SharedRateLimiter(20.0, tmp_path / "rate")
    t0 = time.time()
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(lambda _: lim.acquire(), range(11)))
    assert time.time() - t0 >= 10 / 20.0 * 0.9
    lim.close()


def test_a_429_holds_every_process_back(tmp_path):
    from helpers.ncbi_eutils import SharedRateLimiter
    a = SharedRateLimiter(1000.0, tmp_path / "rate")
    b = SharedRateLimiter(1000.0, tmp_path / "rate")             # "another process"
    a.penalize(0.4)
    t0 = time.time()
    b.acquire()
    assert time.time() - t0 >= 0.35


def test_unusable_state_file_falls_back_to_a_per_process_share(tmp_path, monkeypatch):
    import settings as settings_mod
    from helpers.ncbi_eutils import SharedRateLimiter
    blocker = tmp_path / "file"
    blocker.write_text("x")
    monkeypatch.setattr(settings_mod.settings, "NCBI_PROCESSES", 2)
    lim = SharedRateLimiter(10.0, blocker / "sub" / "rate")      # parent is a file: cannot create
    assert not lim.shared and abs(lim.interval - 0.2) < 1e-9     # 10/2 per process
    lim.acquire()


def test_corrupt_state_does_not_stall(tmp_path):
    import struct
    from helpers.ncbi_eutils import SharedRateLimiter
    path = tmp_path / "rate"
    path.write_bytes(struct.pack("d", time.time() + 10**6))      # far future
    lim = SharedRateLimiter(1000.0, path)
    t0 = time.time()
    lim.acquire()
    assert time.time() - t0 < 1.0
