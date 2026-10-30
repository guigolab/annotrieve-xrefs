"""Runtime settings for the annotrieve-xrefs API."""
from __future__ import annotations

import os
from pathlib import Path


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None or not str(raw).strip():
        return default
    return int(raw)


def _env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None or not str(raw).strip():
        return default
    return float(raw)


class Settings:
    # Published gene_corpus work-dir: {GENE_CORPUS_DIR}/gene_corpus.sqlite
    # and {GENE_CORPUS_DIR}/per_gff/<id>/genes.sqlite
    GENE_CORPUS_DIR: str = os.getenv("GENE_CORPUS_DIR", "")

    # Negative cache_size is KiB (~256 MiB).
    READ_CACHE_SIZE: int = _env_int("READ_CACHE_SIZE", -262144)

    # On-disk annotation GFFs (same host tree as Annotrieve).
    ANNOTATIONS_DIR: str = os.getenv("ANNOTATIONS_DIR", "")
    # Optional harvest --files-root prefix to strip from annotation.gff_path.
    ANNOTATIONS_GFF_PATH_PREFIX: str = os.getenv(
        "ANNOTATIONS_GFF_PATH_PREFIX", ""
    )

    # NCBI E-utilities (key from .env only — never compose literals).
    NCBI_API_KEY: str = os.getenv("NCBI_API_KEY", "")
    NCBI_EMAIL: str = os.getenv("NCBI_EMAIL", "annotrieve-xrefs@crg.eu")
    NCBI_TOOL: str = os.getenv("NCBI_TOOL", "annotrieve-xrefs")
    NCBI_WORKERS: int = _env_int("NCBI_WORKERS", 2)

    # Sequences endpoint.
    SEQUENCES_MAX_WINDOW_BP: int = _env_int("SEQUENCES_MAX_WINDOW_BP", 5_000_000)
    SEQUENCES_COALESCE_GAP_BP: int = _env_int("SEQUENCES_COALESCE_GAP_BP", 10_000)
    SEQUENCES_CACHE_MB: int = _env_int("SEQUENCES_CACHE_MB", 32)
    SEQUENCES_MAX_INFLIGHT: int = _env_int("SEQUENCES_MAX_INFLIGHT", 2)
    SEQUENCES_DEADLINE_S: float = _env_float("SEQUENCES_DEADLINE_S", 800.0)
    GFFREAD_BIN: str = os.getenv("GFFREAD_BIN", "gffread")
    GFFREAD_TIMEOUT_S: float = _env_float("GFFREAD_TIMEOUT_S", 120.0)


settings = Settings()


def gene_corpus_dir() -> Path | None:
    raw = (settings.GENE_CORPUS_DIR or "").strip()
    if not raw:
        return None
    return Path(raw)


def annotations_dir() -> Path | None:
    raw = (settings.ANNOTATIONS_DIR or "").strip()
    if not raw:
        return None
    return Path(raw)
