"""Process-local concurrency gates for expensive API handlers."""
from __future__ import annotations

import threading
from contextlib import contextmanager
from collections.abc import Iterator

# Cap overlapping /hits work per gunicorn worker (RAM back-pressure).
HITS_MAX_INFLIGHT = 16
_HITS_SEMAPHORE = threading.Semaphore(HITS_MAX_INFLIGHT)


class BusyError(Exception):
    """Raised when the hits concurrency gate is saturated."""


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
