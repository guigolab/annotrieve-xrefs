"""Runtime settings for the annotrieve-xrefs API."""
from __future__ import annotations

import os
from pathlib import Path


class Settings:
    # Published gene_corpus work-dir: {GENE_CORPUS_DIR}/gene_corpus.sqlite
    # and {GENE_CORPUS_DIR}/per_gff/<id>/genes.sqlite
    GENE_CORPUS_DIR: str = os.getenv("GENE_CORPUS_DIR", "")

    # Negative cache_size is KiB (~256 MiB).
    READ_CACHE_SIZE: int = int(os.getenv("READ_CACHE_SIZE", "-262144"))


settings = Settings()


def gene_corpus_dir() -> Path | None:
    raw = (settings.GENE_CORPUS_DIR or "").strip()
    if not raw:
        return None
    return Path(raw)
