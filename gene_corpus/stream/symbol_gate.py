"""Name / symbol gating helpers for gene_corpus profiles."""
from __future__ import annotations

import re

# Locus-shaped / ID-echo patterns (Name audit + source map locked gate).
_LOCUS_SHAPED = re.compile(
    r"(?i)^(?:"
    r"LOC\d+"
    r"|.*_LOCUS.*"
    r"|[A-Za-z]+_\d+$"  # e.g. Tco_12345
    r"|gene[-:]?\d+"
    r")$"
)


def is_locus_shaped(value: str) -> bool:
    text = (value or "").strip()
    if not text:
        return True
    return _LOCUS_SHAPED.match(text) is not None


def symbol_ok(
    value: str,
    *,
    locus_tag: str | None = None,
    feature_id: str | None = None,
) -> bool:
    """True when value may be emitted as symbol (not locus/ID echo)."""
    text = (value or "").strip()
    if not text:
        return False
    if locus_tag and text == locus_tag.strip():
        return False
    if feature_id and text == feature_id.strip():
        return False
    if is_locus_shaped(text):
        return False
    return True


def alias_ok(
    value: str,
    *,
    locus_tag: str | None = None,
    feature_id: str | None = None,
) -> bool:
    """
    True when value may be emitted as alias.

    Same gates as symbol_ok, plus reject any whitespace (organism-style
    gb-synonym / old-name strings).
    """
    text = (value or "").strip()
    if not text:
        return False
    if any(ch.isspace() for ch in text):
        return False
    return symbol_ok(text, locus_tag=locus_tag, feature_id=feature_id)


def strip_wrapping_quotes(value: str) -> str:
    """
    Remove one matching pair of wrapping ASCII quotes if present.

    GFF3 says attribute values should not be quoted, but some producers wrap
    Name/gene in ``'...'`` or ``"..."``. One layer only.
    """
    text = value or ""
    if len(text) >= 3 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text


def normalize_symbol(value: str) -> str:
    return strip_wrapping_quotes((value or "").strip()).casefold()
