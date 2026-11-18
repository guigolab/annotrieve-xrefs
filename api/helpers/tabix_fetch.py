"""Thread-local tabix handles and gene-window fetches for GFF files (pysam)."""
from __future__ import annotations

import threading
from collections import OrderedDict
from pathlib import Path

import pysam

from helpers.annotation_files import AnnotationFiles


class FetchError(Exception):
    """A GFF window could not be read. ``code`` is a skip reason / API error code."""

    def __init__(self, code: str, message: str, **extra) -> None:
        super().__init__(message)
        self.code = code
        self.extra = extra


_TABIX_LOCAL = threading.local()
_TABIX_MAX_HANDLES = 32
_TABIX_REGISTRY_LOCK = threading.Lock()
# thread ident -> OrderedDict of path -> TabixFile (for process-wide close)
_TABIX_REGISTRY: dict[int, OrderedDict[str, pysam.TabixFile]] = {}


def open_tabix(gff_path: Path, index_path: Path | None = None):
    """Thread-local LRU of TabixFile handles (pysam is not thread-safe).

    Annotrieve indexes with ``tabix --csi`` (``.csi``), not ``.tbi``. Pass the
    resolved index path so pysam does not default to looking for ``.tbi``.
    """
    cache: OrderedDict[str, pysam.TabixFile] | None = getattr(
        _TABIX_LOCAL, "handles", None
    )
    if cache is None:
        cache = OrderedDict()
        _TABIX_LOCAL.handles = cache
        with _TABIX_REGISTRY_LOCK:
            _TABIX_REGISTRY[threading.get_ident()] = cache
    index = str(index_path) if index_path is not None else None
    key = f"{gff_path}\0{index or ''}"
    handle = cache.get(key)
    if handle is not None:
        cache.move_to_end(key)
        return handle
    handle = (
        pysam.TabixFile(str(gff_path), index=index)
        if index
        else pysam.TabixFile(str(gff_path))
    )
    cache[key] = handle
    while len(cache) > _TABIX_MAX_HANDLES:
        _old_key, old = cache.popitem(last=False)
        try:
            old.close()
        except Exception:
            pass
    return handle


def close_tabix_handles() -> None:
    """Close all thread-local tabix handles (tests / shutdown)."""
    with _TABIX_REGISTRY_LOCK:
        caches = list(_TABIX_REGISTRY.values())
        _TABIX_REGISTRY.clear()
    for cache in caches:
        while cache:
            _key, handle = cache.popitem(last=False)
            try:
                handle.close()
            except Exception:
                pass
    if getattr(_TABIX_LOCAL, "handles", None) is not None:
        _TABIX_LOCAL.handles = None


def fetch_region_lines(
    files: AnnotationFiles,
    *,
    seqid: str,
    start: int,
    end: int,
) -> list[str]:
    # pysam regions are 0-based half-open.
    start0 = max(0, int(start) - 1)
    end0 = max(start0 + 1, int(end))
    try:
        tbx = open_tabix(files.gff_path, files.csi_path)
        return list(tbx.fetch(seqid, start0, end0))
    except ValueError as exc:
        # Unknown contig etc.
        raise FetchError("seqid_not_resolved", f"tabix fetch failed for seqid {seqid!r}: {exc}", seqid=seqid) from exc
    except OSError as exc:
        raise FetchError("annotation_file_not_found", f"failed to read GFF: {exc}", gff_path=str(files.gff_path)) from exc
    except Exception as exc:
        # pysam.SamtoolsError and other htslib failures are bare Exception.
        raise FetchError("annotation_file_not_found", f"failed to read GFF: {exc}", gff_path=str(files.gff_path)) from exc
