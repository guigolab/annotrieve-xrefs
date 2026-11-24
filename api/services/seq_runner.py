"""The single queue runner of ``/downloads``.

Every gunicorn worker starts one runner thread. A file lock (``DOWNLOADS_LOCK_FILE``, on the
container's own disk: all workers share the container) lets only one of them run a download at a
time, so one download at a time uses the NCBI budget and the memory of the extraction. The runner
that holds the lock takes the oldest queued download, runs it, and goes on until none is left.
"""
from __future__ import annotations

import fcntl
import gzip
import json
import logging
import os
import shutil
import tempfile
import threading
import time
from collections import Counter
from pathlib import Path

from helpers.cursor import QueryError
from services.seq_export import run_export
from services.seq_request import parse_request
from services.seq_resolve import SeqError
from services.seq_results import (
    fasta_path,
    oldest_queued,
    purge_expired,
    read_status,
    result_dir,
    write_status,
)
from services.seq_source import Source
from settings import settings

logger = logging.getLogger(__name__)

_thread: threading.Thread | None = None
_stop = threading.Event()


class Cancelled(Exception):
    pass


def _lock_path() -> Path:
    raw = (settings.DOWNLOADS_LOCK_FILE or "").strip()
    return Path(raw) if raw else Path(tempfile.gettempdir()) / "annotrieve-downloads.lock"


def _gzip_fasta(d: Path) -> int:
    """result.fa -> result.fa.gz in blocks (memory flat); returns the compressed size."""
    src, dst = d / "result.fa", d / "result.fa.gz"
    tmp = d / "result.fa.gz.part"
    with open(src, "rb") as fin, gzip.open(tmp, "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, 1024 * 1024)
    os.replace(tmp, dst)
    src.unlink()
    return dst.stat().st_size


def run_download(rid: str) -> None:
    """Run one queued download to the end; its status tells how it went."""
    d = result_dir(rid)
    t0 = time.monotonic()
    try:
        write_status(rid, state="running", started=time.time())
        body = json.loads((d / "request.json").read_text())
        req = parse_request(body)
        src = Source(req)
        try:
            total = (read_status(rid) or {}).get("total")

            def progress(done: int, total_: int | None, counts: Counter) -> None:
                if (d / "cancel").exists():
                    raise Cancelled()
                write_status(rid, state="running", done=done, total=total_ or total,
                             n_ok=counts.get("ok", 0), n_skipped=counts.get("skipped", 0))

            result = run_export(req, src.items(), out_dir=d, total=total, progress=progress)
            data_version = src.data_version
        finally:
            src.close()
        size = _gzip_fasta(d)
        write_status(rid, state="done", n_sequences=result.n_sequences, n_genes=result.n_genes,
                     n_skipped=result.n_skipped, reasons=result.reasons, size_bytes=size,
                     data_version=data_version, done=result.n_genes, error=None)
        logger.info("download %s genes=%d seqs=%d skipped=%d bytes=%d wall_s=%.1f", rid,
                    result.n_genes, result.n_sequences, result.n_skipped, size,
                    time.monotonic() - t0)
    except Cancelled:
        shutil.rmtree(d, ignore_errors=True)
        logger.info("download %s cancelled", rid)
    except (SeqError, QueryError) as exc:
        code = getattr(exc, "code", "download_failed")
        write_status(rid, state="failed", error={
            "code": code, "message": str(exc), "status": getattr(exc, "status", 422)})
        _drop_files(d)
    except Exception as exc:                                        # noqa: BLE001
        logger.exception("download %s failed", rid)
        write_status(rid, state="failed", error={
            "code": "download_failed", "message": str(exc) or exc.__class__.__name__, "status": 500})
        _drop_files(d)


def _drop_files(d: Path) -> None:
    for p in d.iterdir():
        if p.name not in ("status.json", "request.json"):
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            else:
                p.unlink(missing_ok=True)


def run_pending_once(*, stoppable: bool = False) -> int:
    """Take the lock, run every queued download (oldest first). Returns how many were run."""
    path = _lock_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    ran = 0
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return 0                                    # another worker is the runner right now
        try:
            while not (stoppable and _stop.is_set()):
                rid = oldest_queued()
                if rid is None:
                    break
                run_download(rid)
                ran += 1
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
    return ran


def _loop() -> None:
    last_purge = 0.0
    while not _stop.is_set():
        try:
            run_pending_once(stoppable=True)
            if time.monotonic() - last_purge > 60:
                purge_expired()
                last_purge = time.monotonic()
        except Exception:                                           # noqa: BLE001
            logger.exception("downloads runner error")
        _stop.wait(float(settings.DOWNLOADS_POLL_S))


def start_runner() -> None:
    global _thread
    if _thread is not None and _thread.is_alive():
        return
    _stop.clear()
    _thread = threading.Thread(target=_loop, name="downloads-runner", daemon=True)
    _thread.start()


def stop_runner() -> None:
    global _thread
    _stop.set()
    if _thread is not None:
        _thread.join(timeout=5)
    _thread = None
