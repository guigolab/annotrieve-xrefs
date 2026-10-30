"""NCBI E-utilities efetch client with process-wide rate limiting and cache."""
from __future__ import annotations

import logging
import random
import re
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass

import httpx

from helpers.windows import GenomicWindow
from settings import settings

logger = logging.getLogger(__name__)

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

_ACCESSION_RE = re.compile(
    r"^(?:"
    r"[A-Z]{1,2}_\d+\.\d+"
    r"|[A-Z]{1,2}\d{5,}\.\d+"
    r")$"
)

# NCBI: 10 req/s with key, 3 without. Stay under.
_RATE_WITH_KEY = 9.0
_RATE_WITHOUT_KEY = 2.5


class NcbiError(Exception):
    """NCBI fetch failure (maps to API error codes)."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        status: int = 502,
        retry_after: str | None = None,
        **extra,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.status = status
        self.retry_after = retry_after
        self.extra = extra


@dataclass(frozen=True)
class FastaSlice:
    accession: str
    start: int
    stop: int
    sequence: str
    expected_length: int

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def truncated(self) -> bool:
        return self.length != self.expected_length


class _TokenBucket:
    """Simple process-wide token bucket (one token per request)."""

    def __init__(self, rate_per_sec: float) -> None:
        self._interval = 1.0 / rate_per_sec if rate_per_sec > 0 else 0.0
        self._lock = threading.Lock()
        self._next = 0.0

    @property
    def interval(self) -> float:
        return self._interval

    def acquire(self) -> None:
        if self._interval <= 0:
            return
        with self._lock:
            now = time.monotonic()
            if now < self._next:
                delay = self._next - now
                self._next += self._interval
            else:
                delay = 0.0
                self._next = now + self._interval
        if delay > 0:
            time.sleep(delay)


class _WindowCache:
    """Bounded LRU of fetched sequences keyed by (accession, start, stop)."""

    def __init__(self, max_bytes: int) -> None:
        self._max = max(0, int(max_bytes))
        self._lock = threading.Lock()
        self._data: OrderedDict[tuple[str, int, int], str] = OrderedDict()
        self._nbytes = 0

    def get(self, key: tuple[str, int, int]) -> str | None:
        if self._max <= 0:
            return None
        with self._lock:
            seq = self._data.get(key)
            if seq is None:
                return None
            self._data.move_to_end(key)
            return seq

    def put(self, key: tuple[str, int, int], sequence: str) -> None:
        if self._max <= 0:
            return
        size = len(sequence)
        with self._lock:
            old = self._data.pop(key, None)
            if old is not None:
                self._nbytes -= len(old)
            while self._data and self._nbytes + size > self._max:
                _, evicted = self._data.popitem(last=False)
                self._nbytes -= len(evicted)
            if size <= self._max:
                self._data[key] = sequence
                self._nbytes += size

    def clear(self) -> None:
        with self._lock:
            self._data.clear()
            self._nbytes = 0


def _per_process_rate() -> float:
    workers = max(1, int(settings.NCBI_WORKERS))
    key = (settings.NCBI_API_KEY or "").strip()
    total = _RATE_WITH_KEY if key else _RATE_WITHOUT_KEY
    return total / workers


def _parse_fasta_body(text: str) -> str:
    """Return concatenated sequence from a single FASTA record body."""
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise NcbiError(
            "ncbi_error",
            "NCBI response is not FASTA",
            status=502,
        )
    # Reject multi-record or HTML/XML masquerading as 200.
    header_count = sum(1 for line in lines if line.startswith(">"))
    if header_count != 1:
        raise NcbiError(
            "ncbi_error",
            f"NCBI response has {header_count} FASTA records (expected 1)",
            status=502,
        )
    lowered = text[:200].lower()
    if "<html" in lowered or "<?xml" in lowered or "<!doctype" in lowered:
        raise NcbiError(
            "ncbi_error",
            "NCBI returned an HTML/XML error body",
            status=502,
        )
    seq_parts: list[str] = []
    for line in lines[1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            break
        seq_parts.append(s)
    sequence = "".join(seq_parts)
    if not sequence:
        raise NcbiError(
            "ncbi_error",
            "NCBI FASTA record has empty sequence",
            status=502,
        )
    return sequence


class NcbiEutilsClient:
    """Shared httpx client + rate limiter for nuccore slice fetches."""

    def __init__(
        self,
        *,
        rate_per_sec: float | None = None,
        cache_mb: int | None = None,
        transport: httpx.BaseTransport | None = None,
    ) -> None:
        key = (settings.NCBI_API_KEY or "").strip()
        self._api_key = key or None
        self._email = (settings.NCBI_EMAIL or "").strip() or "annotrieve-xrefs@crg.eu"
        self._tool = (settings.NCBI_TOOL or "").strip() or "annotrieve-xrefs"
        rate = rate_per_sec if rate_per_sec is not None else _per_process_rate()
        self._bucket = _TokenBucket(rate)
        mb = settings.SEQUENCES_CACHE_MB if cache_mb is None else cache_mb
        self._cache = _WindowCache(max(0, int(mb)) * 1024 * 1024)
        self._client = httpx.Client(
            timeout=httpx.Timeout(60.0, connect=5.0),
            transport=transport,
            headers={"User-Agent": f"{self._tool}/1.0 ({self._email})"},
        )
        self._closed = False
        logger.info(
            "NCBI efetch client ready: rate=%.2f/s key=%s cache_mb=%s",
            rate,
            "yes" if self._api_key else "no",
            mb,
        )

    @property
    def rate_interval(self) -> float:
        return self._bucket.interval

    @property
    def has_api_key(self) -> bool:
        return self._api_key is not None

    def close(self) -> None:
        if not self._closed:
            self._client.close()
            self._closed = True

    def clear_cache(self) -> None:
        self._cache.clear()

    def fetch_window(
        self,
        window: GenomicWindow,
        *,
        allow_truncation: bool = False,
    ) -> FastaSlice:
        """
        Fetch ``[start, stop]`` (1-based, plus strand) for ``accession``.

        When ``allow_truncation`` is False, a short sequence raises
        ``ncbi_sequence_truncated``. Callers that expect chromosome-end
        clipping pass True and validate the shortfall themselves.
        """
        accession = window.accession.strip()
        if not _ACCESSION_RE.match(accession):
            raise NcbiError(
                "ncbi_error",
                f"invalid INSDC accession {accession!r}",
                status=502,
                accession=accession,
            )
        expected = window.length
        key = (accession, window.start, window.stop)
        cached = self._cache.get(key)
        if cached is not None:
            return FastaSlice(
                accession=accession,
                start=window.start,
                stop=window.stop,
                sequence=cached,
                expected_length=expected,
            )

        params: dict[str, str | int] = {
            "db": "nuccore",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
            "seq_start": window.start,
            "seq_stop": window.stop,
            "strand": 1,
            "tool": self._tool,
            "email": self._email,
        }
        if self._api_key:
            params["api_key"] = self._api_key

        text = self._request_with_retries(params, accession=accession)
        sequence = _parse_fasta_body(text)
        if len(sequence) != expected and not allow_truncation:
            raise NcbiError(
                "ncbi_sequence_truncated",
                (
                    f"NCBI returned {len(sequence)} bp for "
                    f"{accession}:{window.start}-{window.stop} "
                    f"(expected {expected})"
                ),
                status=502,
                accession=accession,
                start=window.start,
                stop=window.stop,
                expected=expected,
                got=len(sequence),
            )
        if len(sequence) == expected:
            self._cache.put(key, sequence)
        return FastaSlice(
            accession=accession,
            start=window.start,
            stop=window.stop,
            sequence=sequence,
            expected_length=expected,
        )

    def _request_with_retries(
        self,
        params: dict[str, str | int],
        *,
        accession: str,
        attempts: int = 3,
    ) -> str:
        last_exc: Exception | None = None
        for attempt in range(attempts):
            self._bucket.acquire()
            try:
                response = self._client.get(EUTILS_EFETCH, params=params)
            except httpx.TimeoutException as exc:
                last_exc = exc
                if attempt + 1 >= attempts:
                    raise NcbiError(
                        "ncbi_timeout",
                        f"NCBI efetch timed out for {accession}",
                        status=504,
                        accession=accession,
                    ) from exc
                time.sleep((2**attempt) * 0.25 + random.uniform(0, 0.1))
                continue
            except httpx.HTTPError as exc:
                last_exc = exc
                if attempt + 1 >= attempts:
                    raise NcbiError(
                        "ncbi_error",
                        f"NCBI efetch failed for {accession}: {exc}",
                        status=502,
                        accession=accession,
                    ) from exc
                time.sleep((2**attempt) * 0.25 + random.uniform(0, 0.1))
                continue

            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After", "1")
                if attempt + 1 >= attempts:
                    raise NcbiError(
                        "ncbi_rate_limited",
                        f"NCBI rate limited for {accession}",
                        status=503,
                        retry_after=str(retry_after),
                        accession=accession,
                    )
                try:
                    delay = float(retry_after)
                except ValueError:
                    delay = (2**attempt) + random.uniform(0, 0.25)
                time.sleep(max(0.1, delay))
                continue

            if response.status_code >= 500:
                if attempt + 1 >= attempts:
                    raise NcbiError(
                        "ncbi_error",
                        f"NCBI server error {response.status_code} for {accession}",
                        status=502,
                        accession=accession,
                        http_status=response.status_code,
                    )
                time.sleep((2**attempt) * 0.5 + random.uniform(0, 0.1))
                continue

            if response.status_code == 404:
                raise NcbiError(
                    "ncbi_sequence_not_found",
                    f"NCBI sequence not found: {accession}",
                    status=502,
                    accession=accession,
                )

            if response.status_code != 200:
                # 4xx other than 429: do not retry.
                body_preview = (response.text or "")[:200]
                raise NcbiError(
                    "ncbi_error",
                    (
                        f"NCBI efetch HTTP {response.status_code} "
                        f"for {accession}: {body_preview}"
                    ),
                    status=502,
                    accession=accession,
                    http_status=response.status_code,
                )

            text = response.text or ""
            # Some NCBI errors arrive as HTTP 200 with an explanatory body.
            if not text.lstrip().startswith(">"):
                preview = text.strip()[:200]
                lowered = preview.lower()
                if "not found" in lowered or "failed" in lowered:
                    raise NcbiError(
                        "ncbi_sequence_not_found",
                        f"NCBI sequence not found: {accession}",
                        status=502,
                        accession=accession,
                    )
                raise NcbiError(
                    "ncbi_error",
                    f"NCBI non-FASTA body for {accession}: {preview}",
                    status=502,
                    accession=accession,
                )
            return text

        raise NcbiError(
            "ncbi_error",
            f"NCBI efetch failed for {accession}",
            status=502,
            accession=accession,
        ) from last_exc


_CLIENT_LOCK = threading.Lock()
_CLIENT: NcbiEutilsClient | None = None


def get_ncbi_client() -> NcbiEutilsClient:
    global _CLIENT
    with _CLIENT_LOCK:
        if _CLIENT is None:
            _CLIENT = NcbiEutilsClient()
        return _CLIENT


def reset_ncbi_client(
    client: NcbiEutilsClient | None = None,
) -> NcbiEutilsClient:
    """Replace the process client (tests). Closes the previous one."""
    global _CLIENT
    with _CLIENT_LOCK:
        if _CLIENT is not None:
            _CLIENT.close()
        _CLIENT = client if client is not None else NcbiEutilsClient()
        return _CLIENT


def close_ncbi_client() -> None:
    global _CLIENT
    with _CLIENT_LOCK:
        if _CLIENT is not None:
            _CLIENT.close()
            _CLIENT = None
