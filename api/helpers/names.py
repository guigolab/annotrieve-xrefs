"""Parse gene-name lists and taxid filters for ``/hits``."""
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

    Empty segments are dropped. Raises if none remain or count exceeds cap.
    """
    pieces: Sequence[str] = (
        (raw_values,) if isinstance(raw_values, str) else raw_values
    )
    names: list[str] = []
    for raw in pieces:
        if raw is None:
            continue
        for part in str(raw).split(","):
            text = part.strip()
            if text:
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


def split_name_term(raw: str) -> tuple[str, str, bool]:
    """
    Return ``(input_raw, normalized_name_form, is_prefix)``.

    ``normalized_name_form`` keeps a trailing ``*`` for prefix terms.
    Raises ``QueryError`` for empty fold or too-short prefix stem.
    """
    text = (raw or "").strip()
    if not text:
        raise QueryError(
            "name tokens must be non-empty",
            code="names_required",
        )
    is_prefix = text.endswith("*")
    stem = text[:-1] if is_prefix else text
    folded = normalize_name_token(stem)
    if not folded:
        raise QueryError(
            "name tokens must be non-empty",
            code="names_required",
        )
    if is_prefix and len(folded) < MIN_PREFIX_LEN:
        raise QueryError(
            f"prefix names require at least {MIN_PREFIX_LEN} characters "
            f"before *",
            code="prefix_too_short",
            min_length=MIN_PREFIX_LEN,
        )
    name_form = f"{folded}*" if is_prefix else folded
    return text, name_form, is_prefix


# Re-export caps for services
__all__ = [
    "MAX_PREFIX_EXPAND",
    "MAX_TAXIDS_FILTER",
    "normalize_name_token",
    "parse_name_list",
    "parse_taxid_list",
    "split_name_term",
]
