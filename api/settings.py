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
    # Corpus work-dir: per_gff/<id>/genes.sqlite and serve/ (published index versions).
    GENE_CORPUS_DIR: str = os.getenv("GENE_CORPUS_DIR", "")

    # Presence index: directory holding the ``current`` symlink to
    # ``v<data_version>/presence.sqlite``. /names and /hits read it (symbol vs
    # alias aware); without it they answer 503.
    PRESENCE_DIR: str = os.getenv("PRESENCE_DIR", "")

    # /loci: base URL of the Annotrieve API (used to build ready-made GFF slice links),
    # export row cap, per-request page cache and concurrent streams per worker.
    ANNOTRIEVE_API_BASE: str = os.getenv(
        "ANNOTRIEVE_API_BASE", "https://genome.crg.es/annotrieve/api/v0")
    # Public base of Annotrieve's file server (whole bgzipped GFFs): {base}/{taxid}/{assembly}/{file}.
    ANNOTRIEVE_FILES_BASE: str = os.getenv(
        "ANNOTRIEVE_FILES_BASE", "https://genome.crg.es/annotrieve/files")
    LOCI_MAX_ROWS: int = _env_int("LOCI_MAX_ROWS", 2_000_000)
    LOCI_CACHE_KIB: int = _env_int("LOCI_CACHE_KIB", 16384)
    LOCI_MAX_INFLIGHT: int = _env_int("LOCI_MAX_INFLIGHT", 3)

    # Negative cache_size is KiB.
    READ_CACHE_SIZE: int = _env_int("READ_CACHE_SIZE", -2048)

    # Search index in RAM. When set, presence.sqlite is copied from the
    # (NFS) version directory into this directory, which must be on tmpfs (for example /dev/shm/...),
    # once per data version and shared by all workers. Empty = read it in place, as before.
    PRESENCE_RAM_DIR: str = os.getenv("PRESENCE_RAM_DIR", "")
    # Largest presence.sqlite (MB) that may be copied into RAM; bigger files are read in place.
    PRESENCE_RAM_MAX_MB: int = _env_int("PRESENCE_RAM_MAX_MB", 450)
    # SQLite page cache (KiB) of a connection that reads the RAM copy: pages are memory-mapped from
    # the tmpfs file, so the private cache can stay small.
    PRESENCE_RAM_CACHE_KIB: int = _env_int("PRESENCE_RAM_CACHE_KIB", 2048)
    # Seconds between checks of the published `current` link (0 = check on every request).
    PRESENCE_RECHECK_S: float = _env_float("PRESENCE_RECHECK_S", 0.0)

    # Threads that run request handlers per worker (framework default 40). Each thread keeps its own
    # read connections, so this also caps the number of open connections.
    THREADPOOL_SIZE: int = _env_int("THREADPOOL_SIZE", 12)

    # /annotations/*: concurrent requests per worker (each opens a multi-GB shard) and the longest
    # a single request may spend reading it (seconds) before it is stopped with a 503.
    ANNOTATIONS_MAX_INFLIGHT: int = _env_int("ANNOTATIONS_MAX_INFLIGHT", 3)
    ANNOTATIONS_QUERY_BUDGET_S: float = _env_float("ANNOTATIONS_QUERY_BUDGET_S", 50.0)

    # Concurrent /hits requests per worker (each can hold tens of MB while it builds its answer).
    HITS_MAX_INFLIGHT: int = _env_int("HITS_MAX_INFLIGHT", 5)

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
    # All API processes on a host share ONE NCBI budget through a small lock file (default: next to
    # the sequences results; set this to put it elsewhere, for example a volume shared by replicas).
    NCBI_RATE_STATE_FILE: str = os.getenv("NCBI_RATE_STATE_FILE", "")
    # Only used when that file cannot be opened: each process then limits itself to
    # rate / NCBI_PROCESSES (set it to the number of gunicorn workers).
    NCBI_PROCESSES: int = _env_int("NCBI_PROCESSES", 1)

    # Downloads (FASTA) endpoint.
    # A download (POST /downloads) holds at most this many genes and this many sequences (a gene
    # gives one sequence unless the transcript rule keeps several).
    SEQUENCES_MAX_GENES: int = _env_int("SEQUENCES_MAX_GENES", 1000)
    SEQUENCES_MAX_SEQUENCES: int = _env_int("SEQUENCES_MAX_SEQUENCES", 1000)
    SEQUENCES_PLAN_MAX_GENES: int = _env_int("SEQUENCES_PLAN_MAX_GENES", 200)
    SEQUENCES_BATCH_GENES: int = _env_int("SEQUENCES_BATCH_GENES", 25)
    # Downloads (requests, status, gzipped FASTA, quota ledger) live here, shared by all workers.
    SEQUENCES_RESULTS_DIR: str = os.getenv("SEQUENCES_RESULTS_DIR", "")
    # How long files stay: after the first complete download, if never downloaded, if failed.
    DOWNLOAD_KEEP_AFTER_FETCH_S: int = _env_int("DOWNLOAD_KEEP_AFTER_FETCH_S", 600)
    DOWNLOAD_KEEP_UNFETCHED_S: int = _env_int("DOWNLOAD_KEEP_UNFETCHED_S", 7200)
    DOWNLOAD_KEEP_FAILED_S: int = _env_int("DOWNLOAD_KEEP_FAILED_S", 3600)
    # Downloads stored at once (any state) and waiting in the queue at once; over either, a new
    # download is refused with 503. Finished ones are never removed to make room.
    SEQUENCES_MAX_RESULTS: int = _env_int("SEQUENCES_MAX_RESULTS", 100)
    DOWNLOADS_MAX_QUEUED: int = _env_int("DOWNLOADS_MAX_QUEUED", 50)
    # Per client IP: downloads waiting or running now, and new downloads in the last 24 hours
    # (like the GFF upload quota in Annotrieve). A request equal to a stored download is free.
    # 0 disables a limit.
    SEQUENCES_IP_MAX_ACTIVE: int = _env_int("SEQUENCES_IP_MAX_ACTIVE", 2)
    DOWNLOADS_PER_DAY: int = _env_int("DOWNLOADS_PER_DAY", 10)
    # One queue for the whole container: a file lock lets one runner at a time work (NCBI budget).
    DOWNLOADS_LOCK_FILE: str = os.getenv("DOWNLOADS_LOCK_FILE", "")
    DOWNLOADS_POLL_S: float = _env_float("DOWNLOADS_POLL_S", 2.0)
    # Proxies in front of the API that append to X-Forwarded-For (Traefik = 1): the client IP is the
    # entry that many places from the right, so a client cannot choose it by sending the header.
    TRUSTED_PROXY_HOPS: int = _env_int("TRUSTED_PROXY_HOPS", 1)
    # Public base of this API (scheme, host and path prefix) used in links, because the app sits
    # behind a path-rewriting proxy and cannot see its public URL.
    PUBLIC_BASE_URL: str = os.getenv("PUBLIC_BASE_URL", "")
    # Skip genes of annotations whose GFF index (.csi) is bigger than this many MiB on disk:
    # opening it needs about as much memory again (0 = no limit).
    SEQUENCES_MAX_INDEX_MB: float = _env_float("SEQUENCES_MAX_INDEX_MB", 256.0)
    SEQUENCES_MAX_WINDOW_BP: int = _env_int("SEQUENCES_MAX_WINDOW_BP", 5_000_000)
    SEQUENCES_COALESCE_GAP_BP: int = _env_int("SEQUENCES_COALESCE_GAP_BP", 10_000)
    SEQUENCES_CACHE_MB: int = _env_int("SEQUENCES_CACHE_MB", 32)
    SEQUENCES_DEADLINE_S: float = _env_float("SEQUENCES_DEADLINE_S", 1800.0)
    # Data-size limit (MB) of each gffread child; it is also the first process the kernel kills when
    # the container runs out of memory. 0 = no limit.
    GFFREAD_MAX_DATA_MB: int = _env_int("GFFREAD_MAX_DATA_MB", 128)
    # Browsers and the nginx cache may reuse an answer this long; after that they ask again with the
    # data-version ETag and get a cheap 304 while the data has not changed. 0 = no Cache-Control.
    HTTP_CACHE_MAX_AGE_S: int = _env_int("HTTP_CACHE_MAX_AGE_S", 300)
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


def presence_ram_dir() -> Path | None:
    raw = (settings.PRESENCE_RAM_DIR or "").strip()
    return Path(raw) if raw else None
