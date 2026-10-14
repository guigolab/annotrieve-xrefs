"""Shared match-mode and CURIE-list parsing for search routes."""
from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

from helpers.cursor import QueryError

MATCH_ANY = "any"
MATCH_ALL = "all"
MATCH_MODES = frozenset({MATCH_ANY, MATCH_ALL})

MatchMode = Literal["any", "all"]


def parse_match(raw: str | None) -> MatchMode:
    """Validate ``match``; default ``any``."""
    text = MATCH_ANY if raw is None else str(raw).strip().lower()
    if text not in MATCH_MODES:
        raise QueryError(
            "match must be any or all",
            code="invalid_match",
        )
    return text  # type: ignore[return-value]


def parse_curie_list(
    raw_values: str | Sequence[str],
    *,
    max_curies: int,
) -> list[str]:
    """Split comma-separated CURIEs, strip, and drop empty segments."""
    pieces: Sequence[str] = (
        (raw_values,) if isinstance(raw_values, str) else raw_values
    )
    curies: list[str] = []
    for raw in pieces:
        if raw is None:
            continue
        for part in str(raw).split(","):
            text = part.strip()
            if text:
                curies.append(text)
    if not curies:
        raise QueryError(
            "at least one curie is required",
            code="curies_required",
        )
    if len(curies) > max_curies:
        raise QueryError(
            f"curies accepts at most {max_curies} values",
            code="too_many_curies",
            max=max_curies,
            received=len(curies),
        )
    return curies
