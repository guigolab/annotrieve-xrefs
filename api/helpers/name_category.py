"""Classify gene-name accessions for ``GET /names`` suggestions."""
from __future__ import annotations

import re

# First-match-wins order matches the former frontend classifyPrefixName.
_NUMBERED = re.compile(r"[_-]\d")
_HOMELOG = re.compile(r"\.[a-z0-9]")
_PSEUDOGENE_TAIL = re.compile(r"[-_]ps")

NAME_CATEGORIES = (
    "pseudogene",
    "numbered",
    "homeolog",
    "predicted",
    "exact",
    "other",
)

CATEGORY_LABELS: dict[str, str] = {
    "pseudogene": "Pseudogene",
    "numbered": "Numbered",
    "homeolog": "Homeolog",
    "predicted": "Predicted",
    "exact": "Exact",
    "other": "Other",
}


def classify_name(name: str) -> str:
    """Return a category id for ``name`` (case-insensitive)."""
    text = (name or "").casefold()
    if ".ps" in text or _PSEUDOGENE_TAIL.search(text):
        return "pseudogene"
    if _NUMBERED.search(text):
        return "numbered"
    if _HOMELOG.search(text):
        return "homeolog"
    if "predicted" in text:
        return "predicted"
    if not any(ch in text for ch in "._-"):
        return "exact"
    return "other"
