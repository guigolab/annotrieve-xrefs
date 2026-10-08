"""GFF strand codes for the public API."""
from __future__ import annotations


def format_strand(value: int) -> str:
    """Map stored strand int (1 / -1 / 0) to ``+`` / ``-`` / ``.``."""
    if int(value) == 1:
        return "+"
    if int(value) == -1:
        return "-"
    return "."
