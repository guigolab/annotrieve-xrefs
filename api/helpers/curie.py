"""CURIE parsing for the annotrieve-xrefs API (sidecar §4)."""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Allow importing gene_corpus from the sibling package when installed beside api/.
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gene_corpus.namespaces.tiers import TIER_A  # noqa: E402

# Public prefix → SQL namespace (slash forms are not legal CURIE prefixes).
PREFIX_TO_NAMESPACE: dict[str, str] = {
    "swiss-prot": "uniprot/swiss-prot",
    "trembl": "uniprot/trembl",
    "imgt": "imgt/gene-db",
}

# Namespaces whose accession is casefolded at harvest (symbol / alias).
_CASEFOLD_NAMESPACES = frozenset({"symbol", "alias"})

_GO_DIGITS_RE = re.compile(r"(?i)^(?:go[:_]?)?(\d{1,7})$")


class CurieError(ValueError):
    """Raised when a CURIE cannot be parsed or the prefix is unknown."""


def normalize_go_accession(raw: str) -> str | None:
    """Return ``GO:`` + 7 digits, or None if *raw* is not a GO id."""
    text = (raw or "").strip()
    if not text:
        return None
    m = _GO_DIGITS_RE.fullmatch(text)
    if m is None:
        return None
    return f"GO:{m.group(1).zfill(7)}"


def resolve_namespace(prefix: str) -> str | None:
    """
    Map a CURIE prefix to a stored namespace.

    Returns None for an unknown prefix.
    """
    key = prefix.strip().lower()
    if not key:
        return None
    if key in PREFIX_TO_NAMESPACE:
        return PREFIX_TO_NAMESPACE[key]
    if key in TIER_A:
        return key
    return None


def normalize_accession(namespace: str, accession: str) -> str:
    """Normalize accession the way harvest stored it (sidecar §4)."""
    text = (accession or "").strip()
    if namespace == "go":
        go = normalize_go_accession(text)
        if go is None:
            return text
        return go
    if namespace in _CASEFOLD_NAMESPACES:
        return text.casefold()
    return text


def parse_curie(raw: str) -> tuple[str, str]:
    """
    Split ``prefix:accession`` on the first colon.

    Returns ``(namespace, accession)``. Raises CurieError when the string
    has no colon or the prefix is empty/unknown.
    """
    text = (raw or "").strip()
    if ":" not in text:
        raise CurieError("CURIE must contain a colon")
    prefix, accession = text.split(":", 1)
    if not prefix.strip():
        raise CurieError("empty CURIE prefix")
    if not accession.strip():
        raise CurieError("empty CURIE accession")
    namespace = resolve_namespace(prefix)
    if namespace is None:
        raise CurieError(f"unknown CURIE prefix: {prefix}")
    return namespace, normalize_accession(namespace, accession)


def format_curie(namespace: str, accession: str) -> str:
    """Public id string for a stored (namespace, accession) pair."""
    for prefix, ns in PREFIX_TO_NAMESPACE.items():
        if ns == namespace:
            return f"{prefix}:{accession}"
    return f"{namespace}:{accession}"
