"""Parse gene-name lists, prefixes, and taxid filters."""
from __future__ import annotations

from collections.abc import Sequence

from helpers.cursor import QueryError
from helpers.pagination import (
    MAX_NAMES_GET,
    MAX_PREFIX_EXPAND,
    MAX_TAXIDS_FILTER,
    MIN_PREFIX_LEN,
)


def normalize_name_token(raw: str) -> str:
    """Strip wrapping quotes and casefold (symbol/alias accession form)."""
    text = (raw or "").strip()
    if len(text) >= 3 and text[0] == text[-1] and text[0] in "'\"":
        text = text[1:-1].strip()
    return text.casefold()


def parse_name_list(
    raw_values: str | Sequence[str],
    *,
    max_names: int = MAX_NAMES_GET,
) -> list[str]:
    """
    Split comma-separated (or sequence) name tokens; preserve order, keep raw.

    Empty segments are dropped. Casefolded duplicates collapse to the first
    spelling before the cap is applied (``too_many_names`` counts uniques).
    """
    pieces: Sequence[str] = (
        (raw_values,) if isinstance(raw_values, str) else raw_values
    )
    names: list[str] = []
    seen: set[str] = set()
    for raw in pieces:
        if raw is None:
            continue
        for part in str(raw).split(","):
            text = part.strip()
            if not text:
                continue
            folded = normalize_name_token(text)
            # Keep raw spelling for terms[].input; fold for uniqueness.
            key = folded or text.casefold()
            if key in seen:
                continue
            seen.add(key)
            names.append(text)
    if not names:
        raise QueryError(
            "at least one name is required",
            code="names_required",
        )
    if len(names) > max_names:
        raise QueryError(
            f"names accepts at most {max_names} values",
            code="too_many_names",
            max=max_names,
            received=len(names),
        )
    return names


def parse_taxid_list(
    raw_values: str | Sequence[str] | Sequence[int] | None,
    *,
    field: str = "taxids",
    max_taxids: int = MAX_TAXIDS_FILTER,
) -> list[int]:
    """
    Parse comma-separated or sequence taxids; positive ints, dedupe, order kept.

    ``None`` or blank → empty list. At most ``max_taxids`` values
    (``too_many_taxids``).
    """
    if raw_values is None:
        return []
    pieces: Sequence = (
        (raw_values,) if isinstance(raw_values, str) else raw_values
    )
    out: list[int] = []
    seen: set[int] = set()
    for raw in pieces:
        if raw is None:
            continue
        if isinstance(raw, int):
            parts = [raw]
        else:
            parts = str(raw).split(",")
        for part in parts:
            if isinstance(part, int):
                taxid = part
                text = str(part)
            else:
                text = part.strip()
                if not text:
                    continue
                try:
                    taxid = int(text)
                except ValueError as exc:
                    raise QueryError(
                        f"{field} values must be positive integers",
                        code="invalid_taxid",
                        field=field,
                        value=text,
                    ) from exc
            if taxid < 1:
                raise QueryError(
                    f"{field} values must be positive integers",
                    code="invalid_taxid",
                    field=field,
                    value=taxid,
                )
            if taxid not in seen:
                seen.add(taxid)
                out.append(taxid)
    if len(out) > max_taxids:
        raise QueryError(
            f"{field} accepts at most {max_taxids} values",
            code="too_many_taxids",
            field=field,
            max=max_taxids,
            received=len(out),
        )
    return out


def parse_exact_name(raw: str) -> tuple[str, str]:
    """
    Return ``(input_raw, normalized_name)`` for an exact gene name.

    Raises ``wildcard_not_supported`` when the token contains ``*``.
    """
    text = (raw or "").strip()
    if not text:
        raise QueryError(
            "name tokens must be non-empty",
            code="names_required",
        )
    if "*" in text:
        raise QueryError(
            "wildcards are not supported on this endpoint; use GET /names?prefix=",
            code="wildcard_not_supported",
        )
    folded = normalize_name_token(text)
    if not folded:
        raise QueryError(
            "name tokens must be non-empty",
            code="names_required",
        )
    return text, folded


def parse_prefix_query(raw: str) -> str:
    """
    Normalize a ``GET /names`` prefix (no wildcard character).

    Strips quotes, casefolds, rejects ``*``, and requires ``MIN_PREFIX_LEN``.
    """
    text = (raw or "").strip()
    if not text:
        raise QueryError(
            "prefix is required",
            code="prefix_too_short",
            min_length=MIN_PREFIX_LEN,
        )
    if "*" in text:
        raise QueryError(
            "prefix must not contain *; pass the stem only",
            code="invalid_prefix",
        )
    folded = normalize_name_token(text)
    if len(folded) < MIN_PREFIX_LEN:
        raise QueryError(
            f"prefix requires at least {MIN_PREFIX_LEN} characters",
            code="prefix_too_short",
            min_length=MIN_PREFIX_LEN,
        )
    return folded


# Re-export caps for services
__all__ = [
    "MAX_PREFIX_EXPAND",
    "MAX_TAXIDS_FILTER",
    "normalize_name_token",
    "parse_exact_name",
    "parse_name_list",
    "parse_prefix_query",
    "parse_taxid_list",
]
