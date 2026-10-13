"""Shared match-mode and CURIE-list parsing for search routes."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from typing import Literal

from helpers.curie import CurieError, format_curie, parse_curie
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


def parse_search_request(
    curies: str | Sequence[str],
    *,
    taxid: int | None,
    match: str | None,
    max_curies: int,
) -> tuple[
    MatchMode,
    list[dict],
    dict[tuple[str, str], list[str]],
    list[str],
    int | None,
]:
    """
    Parse CURIE search inputs.

    Returns ``(match, errors, by_key, display_ids, taxid_i)``.
    """
    match_mode = parse_match(match)
    curie_list = parse_curie_list(curies, max_curies=max_curies)
    if taxid is not None and int(taxid) < 1:
        raise QueryError(
            "taxid must be a positive integer",
            code="invalid_taxid",
        )

    errors: list[dict] = []
    by_key: dict[tuple[str, str], list[str]] = {}
    for raw in curie_list:
        try:
            namespace, accession = parse_curie(raw)
        except CurieError as exc:
            errors.append(
                {"curie": raw, "code": exc.code, "message": str(exc)}
            )
            continue
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )
    display_ids = sorted(
        {format_curie(ns, acc) for ns, acc in by_key}
    )
    taxid_i = int(taxid) if taxid is not None else None
    return match_mode, errors, by_key, display_ids, taxid_i


def taxid_has_annotations(conn: sqlite3.Connection, taxid: int) -> bool:
    """True when any annotation lineage row carries *taxid*."""
    row = conn.execute(
        """
        SELECT 1 FROM annotation_lineage WHERE taxid = ? LIMIT 1
        """,
        (int(taxid),),
    ).fetchone()
    return row is not None
