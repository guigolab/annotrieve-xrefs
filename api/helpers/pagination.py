"""Shared pagination limits and OpenAPI descriptions."""
from __future__ import annotations

import base64
import binascii
import json

DEFAULT_LIMIT = 100
MAX_LIMIT = 200
MAX_OFFSET = 10_000
MAX_CURIES_GET = 20
MAX_ACCESSIONS_FILTER = MAX_CURIES_GET
MIN_Q_LEN = 3

# Gene-name /hits caps (GET stays cacheable via nginx; POST is fallback)
MAX_NAMES_GET = 200
MAX_NAMES_POST = 200
MAX_PREFIX_EXPAND = 500
MIN_PREFIX_LEN = 3
MAX_TAXIDS_FILTER = 50

# POST /sequences (bulk FASTA export)
MAX_SEQUENCES_POST = 200
MAX_SEQUENCES_BUFFER_BP = 10_000
MAX_SEQUENCES_ID_CHARS = 256


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


def encode_offset_cursor(offset: int) -> str:
    """Opaque cursor for ``/annotations/*`` lists: the position of the next page."""
    raw = json.dumps({"o": int(offset)}, separators=(",", ":")).encode()
    return base64.urlsafe_b64encode(raw).decode().rstrip("=")


def resolve_offset(cursor: str | None, offset: int) -> int:
    """The offset to read from: ``cursor`` (from a previous ``next_cursor``) or ``offset``, not both."""
    from helpers.cursor import QueryError

    if cursor is None or cursor == "":
        return offset
    if offset:
        raise QueryError("use either cursor or offset, not both", code="invalid_combination",
                         cursor=cursor, offset=offset)
    try:
        data = json.loads(base64.urlsafe_b64decode(cursor + "=" * (-len(cursor) % 4)))
        value = int(data["o"])
    except (ValueError, KeyError, TypeError, binascii.Error):
        raise QueryError("invalid cursor", code="invalid_cursor") from None
    if value < 0 or value > MAX_OFFSET:
        raise QueryError(f"cursor is outside the pageable range (at most {MAX_OFFSET} rows deep)",
                         code="invalid_cursor")
    return value


def page_envelope(
    *,
    total: int,
    limit: int,
    offset: int,
    results: list,
    **extras,
) -> dict:
    """Standard list envelope: total / limit / offset / next_cursor / results (+ extras)."""
    end = int(offset) + len(results)
    out: dict = {
        "total": int(total),
        "limit": int(limit),
        "offset": int(offset),
        "next_cursor": encode_offset_cursor(end) if end < int(total) and end <= MAX_OFFSET else None,
        "results": results,
    }
    out.update(extras)
    return out
