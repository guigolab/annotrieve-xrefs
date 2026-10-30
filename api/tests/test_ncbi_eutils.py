"""NCBI efetch client tests (mocked transport)."""
from __future__ import annotations

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
