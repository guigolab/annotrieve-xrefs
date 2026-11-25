"""Process-local concurrency gates for expensive API handlers."""
from __future__ import annotations

import threading
from contextlib import contextmanager
from collections.abc import Iterator

from settings import settings

# Cap overlapping annotation-search work per gunicorn worker (RAM back-pressure).
HITS_MAX_INFLIGHT = max(1, int(settings.HITS_MAX_INFLIGHT))
_HITS_SEMAPHORE = threading.Semaphore(HITS_MAX_INFLIGHT)


_ANNOTATIONS_SEMAPHORE = threading.Semaphore(max(1, int(settings.ANNOTATIONS_MAX_INFLIGHT)))
_LOCI_SEMAPHORE = threading.Semaphore(max(1, int(settings.LOCI_MAX_INFLIGHT)))
_PLAN_SEMAPHORE = threading.Semaphore(max(1, int(settings.SEQUENCES_PLAN_MAX_INFLIGHT)))


class BusyError(Exception):
    """Raised when a concurrency gate is saturated."""


@contextmanager
def hits_slot(*, blocking: bool = False) -> Iterator[None]:
    """
    Acquire a hits in-flight slot.

    Non-blocking by default: raises ``BusyError`` if all slots are taken.
    """
    if not _HITS_SEMAPHORE.acquire(blocking=blocking):
        raise BusyError("too many concurrent hits requests")
    try:
        yield
    finally:
        _HITS_SEMAPHORE.release()


def acquire_loci_slot() -> None:
    """Take a /loci slot (non-blocking); the caller must call ``release_loci_slot``."""
    if not _LOCI_SEMAPHORE.acquire(blocking=False):
        raise BusyError("too many concurrent loci requests")


def release_loci_slot() -> None:
    _LOCI_SEMAPHORE.release()


def annotations_slot() -> Iterator[None]:
    """
    FastAPI dependency for the /annotations/* routes: at most N run at once per worker.

    When all slots are taken the request is refused at once with 503 and ``Retry-After``.
    """
    if not _ANNOTATIONS_SEMAPHORE.acquire(blocking=False):
        from helpers.errors import error_detail
        from fastapi import HTTPException

        raise HTTPException(
            status_code=503,
            detail=error_detail("annotations_busy",
                                "too many concurrent annotation requests; retry shortly"),
            headers={"Retry-After": "1"},
        )
    try:
        yield
    finally:
        _ANNOTATIONS_SEMAPHORE.release()


@contextmanager
def plan_slot() -> Iterator[None]:
    """Gate for the routes that read GFF slices gene by gene (``/transcripts``, the sequence count of
    ``POST /downloads``): at most N at once per worker, the rest are refused at once with 503."""
    if not _PLAN_SEMAPHORE.acquire(blocking=False):
        from helpers.errors import error_detail
        from fastapi import HTTPException

        raise HTTPException(
            status_code=503,
            detail=error_detail("plan_busy", "too many concurrent transcript requests; retry shortly"),
            headers={"Retry-After": "2"},
        )
    try:
        yield
    finally:
        _PLAN_SEMAPHORE.release()


class GuardedBody:
    """Iterator over a streaming body that runs ``cleanup`` exactly once.

    A generator that was never started does not run its ``finally`` when it is closed or collected,
    so a client that goes away before the first chunk would keep a concurrency slot and open
    connections for ever. This wrapper also cleans up when it is closed or garbage collected.
    """

    def __init__(self, body: Iterator, cleanup) -> None:
        self._body = body
        self._cleanup = cleanup
        self._lock = threading.Lock()
        self._done = False

    def __iter__(self) -> "GuardedBody":
        return self

    def __next__(self):
        try:
            return next(self._body)
        except BaseException:
            self.close()
            raise

    def close(self) -> None:
        with self._lock:
            if self._done:
                return
            self._done = True
        try:
            close = getattr(self._body, "close", None)
            if close is not None:
                close()
        finally:
            self._cleanup()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:                                        # noqa: BLE001 - interpreter exit
            pass
