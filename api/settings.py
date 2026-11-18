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

    # Presence index (PLAN-build-tables.md): directory holding the ``current`` symlink to
    # ``v<data_version>/presence.sqlite``. When set, /names and /hits read it (symbol vs
    # alias aware) instead of gene_corpus.sqlite xref_file. Empty = legacy behaviour.
    PRESENCE_DIR: str = os.getenv("PRESENCE_DIR", "")

    # /loci: base URL of the Annotrieve API (used to build ready-made GFF slice links),
    # export row cap, per-request page cache and concurrent streams per worker.
    ANNOTRIEVE_API_BASE: str = os.getenv(
        "ANNOTRIEVE_API_BASE", "https://genome.crg.es/annotrieve/api/v0")
    LOCI_MAX_ROWS: int = _env_int("LOCI_MAX_ROWS", 2_000_000)
    LOCI_CACHE_KIB: int = _env_int("LOCI_CACHE_KIB", 16384)
    LOCI_MAX_INFLIGHT: int = _env_int("LOCI_MAX_INFLIGHT", 4)

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
    # Number of API processes (gunicorn -w): each has its own token bucket, so the NCBI budget is
    # divided by this to keep the whole container within it.
    NCBI_PROCESSES: int = _env_int("NCBI_PROCESSES", 1)

    # Sequences endpoint.
    # Public base of this API (scheme, host and path prefix) used in job/result links, because
    # the app sits behind a path-rewriting proxy and cannot see its public URL.
    PUBLIC_BASE_URL: str = os.getenv("PUBLIC_BASE_URL", "")
    SEQUENCES_MAX_GENES: int = _env_int("SEQUENCES_MAX_GENES", 2000)
    SEQUENCES_SYNC_MAX_GENES: int = _env_int("SEQUENCES_SYNC_MAX_GENES", 50)
    SEQUENCES_PLAN_MAX_GENES: int = _env_int("SEQUENCES_PLAN_MAX_GENES", 200)
    SEQUENCES_BATCH_GENES: int = _env_int("SEQUENCES_BATCH_GENES", 25)
    SEQUENCES_RESULTS_DIR: str = os.getenv("SEQUENCES_RESULTS_DIR", "")
    SEQUENCES_RESULT_TTL_S: int = _env_int("SEQUENCES_RESULT_TTL_S", 86400)
    # Skip genes of annotations whose GFF index (.csi) is bigger than this many MiB on disk:
    # opening it needs about as much memory again (0 = no limit).
    SEQUENCES_MAX_INDEX_MB: float = _env_float("SEQUENCES_MAX_INDEX_MB", 0.0)
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


def presence_dir() -> Path | None:
    raw = (settings.PRESENCE_DIR or "").strip()
    return Path(raw) if raw else None
