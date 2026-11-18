"""NCBI E-utilities efetch client with process-wide rate limiting and cache."""
from __future__ import annotations

import fcntl
import logging
import os
import random
import struct
import re
import tempfile
import threading
import time
from pathlib import Path
from collections import OrderedDict
from dataclasses import dataclass

import httpx

from helpers.windows import GenomicWindow
from settings import settings

logger = logging.getLogger(__name__)

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

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
    """Rate limiter local to this process (one token per request).

    Used for explicit rates (tests) and as the fallback when the shared state file cannot be used.
    """

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

    def penalize(self, seconds: float) -> None:
        """No request before ``seconds`` from now (after NCBI answered 429)."""
        with self._lock:
            self._next = max(self._next, time.monotonic() + max(0.0, seconds))


class SharedRateLimiter:
    """One NCBI budget for every API process on this host.

    The time of the next free request slot lives in a small file guarded by ``flock``. A caller
    reserves the next slot under the lock (a few microseconds), releases it and sleeps until the
    slot, so concurrent processes and threads are spaced ``1/rate`` apart without any of them
    holding the lock while waiting. A 429 from NCBI pushes the slot of everyone forward.
    """

    MAX_AHEAD_S = 120.0                     # a "next" further away than this is a corrupt file

    def __init__(self, rate_per_sec: float, state_file: str | Path) -> None:
        self._interval = 1.0 / rate_per_sec if rate_per_sec > 0 else 0.0
        self._path = Path(state_file)
        self._thread_lock = threading.Lock()
        self._fd: int | None = None
        self._fallback: _TokenBucket | None = None
        try:
            self._path.parent.mkdir(parents=True, exist_ok=True)
            self._fd = os.open(self._path, os.O_RDWR | os.O_CREAT, 0o666)
        except OSError as exc:
            logger.warning("NCBI rate state %s unusable (%s); limiting per process", self._path, exc)
            procs = max(1, int(settings.NCBI_PROCESSES))
            self._fallback = _TokenBucket(rate_per_sec / procs)

    @property
    def interval(self) -> float:
        return self._interval if self._fallback is None else self._fallback.interval

    @property
    def shared(self) -> bool:
        return self._fallback is None

    def _reserve(self, delta: float, *, push_only: bool = False) -> float:
        """Take the next slot (``delta`` later than the previous one); return its wall time."""
        with self._thread_lock:
            fcntl.flock(self._fd, fcntl.LOCK_EX)
            try:
                raw = os.pread(self._fd, 8, 0)
                nxt = struct.unpack("d", raw)[0] if len(raw) == 8 else 0.0
                now = time.time()
                if not (0.0 <= nxt <= now + self.MAX_AHEAD_S):
                    nxt = 0.0
                if push_only:
                    slot = max(nxt, now + delta)
                    new_next = slot
                else:
                    slot = max(nxt, now)
                    new_next = slot + delta
                os.pwrite(self._fd, struct.pack("d", new_next), 0)
                return slot
            finally:
                fcntl.flock(self._fd, fcntl.LOCK_UN)

    def acquire(self) -> None:
        if self._fallback is not None:
            self._fallback.acquire()
            return
        if self._interval <= 0:
            return
        wait = self._reserve(self._interval) - time.time()
        if wait > 0:
            time.sleep(wait)

    def penalize(self, seconds: float) -> None:
        if self._fallback is not None:
            self._fallback.penalize(seconds)
        elif self._interval > 0:
            self._reserve(max(0.0, seconds), push_only=True)

    def close(self) -> None:
        if self._fd is not None:
            try:
                os.close(self._fd)
            except OSError:
                pass
            self._fd = None


def default_state_file() -> Path:
    """Where the shared budget lives: next to the results, so replicas that share that volume share
    the budget too (``NCBI_RATE_STATE_FILE`` overrides)."""
    if settings.NCBI_RATE_STATE_FILE:
        return Path(settings.NCBI_RATE_STATE_FILE)
    root = Path(settings.SEQUENCES_RESULTS_DIR) if settings.SEQUENCES_RESULTS_DIR else (
        Path(tempfile.gettempdir()) / "annotrieve-sequences")
    return root / ".ncbi-rate"


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


def _default_rate() -> float:
    """The whole host's NCBI budget (requests per second); every process shares it through
    :class:`SharedRateLimiter`. ``NCBI_WORKERS`` only sets thread parallelism."""
    key = (settings.NCBI_API_KEY or "").strip()
    return _RATE_WITH_KEY if key else _RATE_WITHOUT_KEY


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
        state_file: str | Path | None = None,
    ) -> None:
        key = (settings.NCBI_API_KEY or "").strip()
        self._api_key = key or None
        self._email = (settings.NCBI_EMAIL or "").strip() or "annotrieve-xrefs@crg.eu"
        self._tool = (settings.NCBI_TOOL or "").strip() or "annotrieve-xrefs"
        rate = rate_per_sec if rate_per_sec is not None else _default_rate()
        # An explicit rate is a local limit (tests, scripts); the service uses the shared budget.
        if rate_per_sec is None or state_file is not None:
            self._bucket = SharedRateLimiter(rate, state_file or default_state_file())
        else:
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
            "NCBI efetch client ready: rate=%.2f/s %s key=%s cache_mb=%s ncbi_workers=%d",
            rate,
            "shared across processes" if getattr(self._bucket, "shared", False)
            else "this process only",
            "yes" if self._api_key else "no",
            mb,
            max(1, int(settings.NCBI_WORKERS)),
        )

    @property
    def rate_interval(self) -> float:
        return self._bucket.interval

    @property
    def has_api_key(self) -> bool:
        return self._api_key is not None

    def close(self) -> None:
        if not self._closed:
            if isinstance(self._bucket, SharedRateLimiter):
                self._bucket.close()
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
        if not accession:
            raise NcbiError(
                "ncbi_error",
                "empty accession",
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
                delay = min(max(0.1, delay), 30.0)
                # NCBI says we are too fast: hold every process back, not just this request.
                self._bucket.penalize(delay)
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
                # NCBI sometimes spaces its letters out ("F a i l e d ..."): compare without blanks.
                compact = re.sub(r"[\s\x00]+", "", preview).lower()
                if "notfound" in compact or "failed" in compact:
                    raise NcbiError(
                        "ncbi_sequence_not_found",
                        f"NCBI has no sequence for {accession}",
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
