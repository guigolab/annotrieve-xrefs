"""Process-local concurrency gates for expensive API handlers."""
from __future__ import annotations

import threading
from contextlib import contextmanager
from collections.abc import Iterator

from settings import settings

# Cap overlapping annotation-search work per gunicorn worker (RAM back-pressure).
HITS_MAX_INFLIGHT = max(1, int(settings.HITS_MAX_INFLIGHT))
_HITS_SEMAPHORE = threading.Semaphore(HITS_MAX_INFLIGHT)

_SEQUENCES_SEMAPHORE = threading.Semaphore(
    max(1, int(settings.SEQUENCES_MAX_INFLIGHT))
)


_ANNOTATIONS_SEMAPHORE = threading.Semaphore(max(1, int(settings.ANNOTATIONS_MAX_INFLIGHT)))
_LOCI_SEMAPHORE = threading.Semaphore(max(1, int(settings.LOCI_MAX_INFLIGHT)))


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


@contextmanager
def sequences_slot(*, blocking: bool = False) -> Iterator[None]:
    """
    Acquire a sequences in-flight slot.

    Non-blocking by default: raises ``BusyError`` if all slots are taken.
    """
    if not _SEQUENCES_SEMAPHORE.acquire(blocking=blocking):
        raise BusyError("too many concurrent sequences requests")
    try:
        yield
    finally:
        _SEQUENCES_SEMAPHORE.release()


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
