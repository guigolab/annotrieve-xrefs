"""Tabix reads of GFF windows (pysam), with handles scoped to one request.

Opening a ``TabixFile`` loads the whole ``.csi`` index into memory (tens of MiB for a large
genome, about 35 ms on NFS). A process- or thread-wide cache of open files therefore grows with
the number of distinct annotations a worker has seen, which is unbounded in practice (a gene
name can occur in thousands of annotations). So a :class:`TabixCache` belongs to one request
(one ``GeneResolver``), keeps at most ``max_handles`` files open (default 1: genes of one
annotation arrive together) and releases them in :meth:`close`.
"""
from __future__ import annotations

import ctypes
from collections import OrderedDict
from pathlib import Path

import pysam

from helpers.annotation_files import AnnotationFiles


def release_memory() -> None:
    """Hand freed heap back to the OS (glibc). Opening a tabix index allocates tens to hundreds of
    MiB for fragmented assemblies; without this the process keeps that memory after the file is
    closed (measured: 170 MiB retained vs 66 MiB after trimming). No-op where unavailable."""
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (OSError, AttributeError):
        pass


class FetchError(Exception):
    """A GFF window could not be read. ``code`` is a skip reason / API error code."""

    def __init__(self, code: str, message: str, **extra) -> None:
        super().__init__(message)
        self.code = code
        self.extra = extra


class TabixCache:
    def __init__(self, max_handles: int = 1) -> None:
        self.max_handles = max(1, int(max_handles))
        self._handles: OrderedDict[str, pysam.TabixFile] = OrderedDict()

    def __len__(self) -> int:
        return len(self._handles)

    def _open(self, gff_path: Path, index_path: Path | None) -> pysam.TabixFile:
        index = str(index_path) if index_path is not None else None
        key = f"{gff_path}\0{index or ''}"
        handle = self._handles.get(key)
        if handle is not None:
            self._handles.move_to_end(key)
            return handle
        handle = (pysam.TabixFile(str(gff_path), index=index) if index
                  else pysam.TabixFile(str(gff_path)))
        self._handles[key] = handle
        while len(self._handles) > self.max_handles:
            _key, old = self._handles.popitem(last=False)
            try:
                old.close()
            except Exception:                                   # noqa: BLE001
                pass
            release_memory()
        return handle

    def fetch(self, files: AnnotationFiles, *, seqid: str, start: int, end: int) -> list[str]:
        """GFF lines overlapping ``seqid:start-end`` (1-based, inclusive)."""
        # pysam regions are 0-based half-open.
        start0 = max(0, int(start) - 1)
        end0 = max(start0 + 1, int(end))
        try:
            return list(self._open(files.gff_path, files.csi_path).fetch(seqid, start0, end0))
        except ValueError as exc:                               # unknown contig
            raise FetchError("seqid_not_resolved",
                             f"tabix fetch failed for seqid {seqid!r}: {exc}", seqid=seqid) from exc
        except Exception as exc:                                # OSError and htslib errors
            raise FetchError("annotation_file_not_found", f"failed to read GFF: {exc}",
                             gff_path=str(files.gff_path)) from exc

    def close(self) -> None:
        had = bool(self._handles)
        while self._handles:
            _key, handle = self._handles.popitem(last=False)
            try:
                handle.close()
            except Exception:                                   # noqa: BLE001
                pass
        if had:
            release_memory()
