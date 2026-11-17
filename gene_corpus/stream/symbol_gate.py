"""Name / symbol gating helpers for gene_corpus profiles."""
from __future__ import annotations

import re

# Locus-shaped / ID-echo patterns (Name audit + source map locked gate).
#
# Two tiers (decision of 2026-11): symbols (Name, gene, standard_name) are only rejected
# for the *core* patterns, so copy-suffixed symbols such as ``Polq_1``, ``Zan_0`` or
# ``selenot_1`` (homeologs, paralog copies) are kept. Aliases keep the strict tier, which
# also rejects ``letters_digits`` shapes, so old_locus_tag-style values do not flood the
# index. A locus tag used as a symbol is still caught when it equals the gene's own
# ``locus_tag`` or ``ID`` (see ``symbol_ok``).
_CORE_LOCUS_SHAPED = re.compile(
    r"(?i)^(?:"
    r"LOC\d+"
    r"|.*_LOCUS.*"
    r"|gene[-:]?\d+"
    r")$"
)
_ALIAS_LOCUS_SHAPED = re.compile(r"(?i)^[A-Za-z]+_\d+$")  # e.g. Tco_12345


def is_locus_shaped_core(value: str) -> bool:
    """Shapes rejected for symbols and aliases alike."""
    text = (value or "").strip()
    if not text:
        return True
    return _CORE_LOCUS_SHAPED.match(text) is not None


def is_locus_shaped(value: str) -> bool:
    """Strict tier (used for aliases): core shapes plus ``letters_digits``."""
    if is_locus_shaped_core(value):
        return True
    return _ALIAS_LOCUS_SHAPED.match((value or "").strip()) is not None


def _echo_or_empty(
    text: str, locus_tag: str | None, feature_id: str | None
) -> bool:
    if not text:
        return True
    if locus_tag and text == locus_tag.strip():
        return True
    if feature_id and text == feature_id.strip():
        return True
    return False


def symbol_ok(
    value: str,
    *,
    locus_tag: str | None = None,
    feature_id: str | None = None,
) -> bool:
    """True when value may be emitted as symbol (not locus/ID echo)."""
    text = (value or "").strip()
    if _echo_or_empty(text, locus_tag, feature_id):
        return False
    return not is_locus_shaped_core(text)


def alias_ok(
    value: str,
    *,
    locus_tag: str | None = None,
    feature_id: str | None = None,
) -> bool:
    """
    True when value may be emitted as alias.

    Strict tier: same echo and core-shape gates as symbol_ok, the extra
    ``letters_digits`` shape, and no whitespace (organism-style gb-synonym /
    old-name strings).
    """
    text = (value or "").strip()
    if _echo_or_empty(text, locus_tag, feature_id):
        return False
    if any(ch.isspace() for ch in text):
        return False
    return not is_locus_shaped(text)


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
