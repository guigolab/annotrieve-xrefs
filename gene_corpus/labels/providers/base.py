"""Label provider protocol for ontology / controlled-vocab loaders."""
from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True)
class LabelTerm:
    """One controlled-vocab term ready to intersect with xref_meta."""

    accession: str
    label: str
    aliases: tuple[str, ...] = ()
    alt_ids: tuple[str, ...] = ()


class LabelProvider(Protocol):
    """Download and stream terms for one corpus namespace."""

    namespace: str

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        """Return a local source file path (cached unless *force*)."""
        ...

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from *path* without loading the whole file into RAM."""
        ...
