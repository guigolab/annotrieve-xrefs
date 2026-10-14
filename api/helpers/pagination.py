"""Shared pagination limits and OpenAPI descriptions."""
from __future__ import annotations

DEFAULT_LIMIT = 100
MAX_LIMIT = 200
MAX_OFFSET = 10_000
MAX_CURIES_GET = 20
MAX_ACCESSIONS_FILTER = MAX_CURIES_GET
MIN_Q_LEN = 3

# Gene-name /hits caps
MAX_NAMES_GET = 50
MAX_NAMES_POST = 200
MAX_PREFIX_EXPAND = 500
MIN_PREFIX_LEN = 3
MAX_TAXIDS_FILTER = 50


def clamp_limit(limit: int) -> int:
    return max(1, min(int(limit), MAX_LIMIT))


def clamp_offset(offset: int) -> int:
    """Return a non-negative offset; raise ValueError if above MAX_OFFSET."""
    value = int(offset)
    if value < 0:
        raise ValueError("offset must be >= 0")
    if value > MAX_OFFSET:
        raise ValueError(f"offset must be <= {MAX_OFFSET}")
    return value


def page_envelope(
    *,
    total: int,
    limit: int,
    offset: int,
    results: list,
    **extras,
) -> dict:
    """Standard list envelope: total / limit / offset / results (+ extras)."""
    out: dict = {
        "total": int(total),
        "limit": int(limit),
        "offset": int(offset),
        "results": results,
    }
    out.update(extras)
    return out
