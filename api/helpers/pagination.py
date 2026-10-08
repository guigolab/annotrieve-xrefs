"""Shared pagination limits and OpenAPI descriptions."""
from __future__ import annotations

DEFAULT_LIMIT = 100
MAX_LIMIT = 200
MAX_CURIES_GET = 20
MAX_CURIES_POST = 100
MAX_ACCESSIONS_FILTER = MAX_CURIES_GET

PAGE_TOKEN_DESC = (
    "Opaque page token from a previous response next/previous field "
    "(mutually exclusive with the other)"
)


def clamp_limit(limit: int) -> int:
    return max(1, min(int(limit), MAX_LIMIT))
