"""The single queue runner of ``/downloads``.

Every gunicorn worker starts one runner thread. A file lock (``DOWNLOADS_LOCK_FILE``, on the
container's own disk: all workers share the container) lets only one of them run a download at a
time, so one download at a time uses the NCBI budget and the memory of the extraction. The runner
that holds the lock takes the oldest queued download, runs it, and goes on until none is left.

Stopping and crashing. A runner that is asked to stop (worker shutdown) puts the running download
back in the queue at its old place. A worker that dies without that chance (kill, out of memory)
leaves a ``running`` status behind; the next runner to take the lock sees it as an orphan, because
nobody else can be running a download then, and queues it again at once instead of waiting for it
to look stale. A download that was orphaned ``MAX_ATTEMPTS`` times is failed (``download_lost``):
it probably kills its worker.

A failure that cannot be recorded (the status file cannot be written: full disk) is counted per
download. The download is held back with a growing delay, and after ``MAX_STRIKES`` it is removed,
so one broken download never keeps the runner busy in a two-second loop.
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
    remove,
    result_dir,
    running_ids,
    write_status,
)
from services.seq_source import Source
from settings import settings

logger = logging.getLogger(__name__)

_thread: threading.Thread | None = None
_stop = threading.Event()


MAX_ATTEMPTS = 3            # runs that ended with the worker dying before a download is given up
MAX_STRIKES = 5             # times recording the outcome of a download failed before it is removed
HOLD_MAX_S = 300.0          # longest wait before a download with a strike is tried again
STOP_WAIT_S = 20.0          # how long stop_runner waits for the running download to step aside

_strikes: dict[str, int] = {}
_held_until: dict[str, float] = {}      # id -> time.monotonic() before which it is not taken


class Cancelled(Exception):
    pass


class Interrupted(Exception):
    """The worker is stopping: the download goes back to the queue and is done later."""


def _lock_path() -> Path:
    raw = (settings.DOWNLOADS_LOCK_FILE or "").strip()
    return Path(raw) if raw else Path(tempfile.gettempdir()) / "annotrieve-downloads.lock"


def _gzip_fasta(d: Path) -> int:
    """result.fa -> result.fa.gz in blocks (memory flat); returns the compressed size."""
    src, dst = d / "result.fa", d / "result.fa.gz"
    tmp = d / "result.fa.gz.part"
    # No name and no time in the gzip header: `gunzip -N` and desktop unzip tools would otherwise
    # name the output after the temp file ("result.fa.gz.part") instead of after the download.
    with open(src, "rb") as fin, open(tmp, "wb") as raw, \
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as fout:
        shutil.copyfileobj(fin, fout, 1024 * 1024)
    os.replace(tmp, dst)
    src.unlink()
    return dst.stat().st_size


def _record_failure(rid: str, error: dict) -> None:
    """Write the ``failed`` status. If the disk refuses, free the download's files and try once more."""
    d = result_dir(rid)
    try:
        write_status(rid, state="failed", error=error)
    except OSError:
        logger.exception("download %s: cannot write the failed status; dropping its files", rid)
        _drop_files(d)
        write_status(rid, state="failed", error=error)
    _drop_files(d)


def run_download(rid: str, *, interruptible: bool = False) -> None:
    """Run one queued download to the end; its status tells how it went.

    With ``interruptible`` a stop request (:func:`stop_runner`) ends the run at the next progress
    report and puts the download back in the queue.
    """
    d = result_dir(rid)
    t0 = time.monotonic()
    try:
        attempts = int((read_status(rid, _raw=True) or {}).get("attempts") or 0)
        write_status(rid, state="running", started=time.time(), attempts=attempts + 1)
        body = json.loads((d / "request.json").read_text())
        req = parse_request(body)
        src = Source(req)
        try:
            total = (read_status(rid) or {}).get("total")

            def progress(done: int, total_: int | None, counts: Counter) -> None:
                if (d / "cancel").exists():
                    raise Cancelled()
                if interruptible and _stop.is_set():
                    raise Interrupted()
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
    except Interrupted:
        # Not the download's fault: give the attempt back and keep its place in the queue.
        _drop_files(d)
        write_status(rid, state="queued", attempts=attempts, done=0, error=None)
        logger.info("download %s put back in the queue: the worker is stopping", rid)
    except (SeqError, QueryError) as exc:
        code = getattr(exc, "code", "download_failed")
        _record_failure(rid, {"code": code, "message": str(exc),
                              "status": getattr(exc, "status", 422)})
    except Exception as exc:                                        # noqa: BLE001
        logger.exception("download %s failed", rid)
        # The cause stays in the log: the status is shown to whoever holds the download id.
        _record_failure(rid, {
            "code": "download_failed", "status": 500,
            "message": "the download failed unexpectedly; ask for it again, or narrow the request"})


def _drop_files(d: Path) -> None:
    for p in d.iterdir():
        if p.name not in ("status.json", "request.json"):
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            else:
                p.unlink(missing_ok=True)


def recover_orphans() -> int:
    """Queue again the downloads a stopped worker left ``running``. Call it holding the lock.

    A download that was cancelled meanwhile is removed; one that already used ``MAX_ATTEMPTS``
    runs is failed. Returns how many were queued again.
    """
    queued = 0
    for rid in running_ids():
        try:
            d = result_dir(rid)
            st = read_status(rid, _raw=True) or {}
            if (d / "cancel").exists():
                shutil.rmtree(d, ignore_errors=True)
            elif int(st.get("attempts") or 0) >= MAX_ATTEMPTS:
                logger.error("download %s was running when its worker stopped %d times: failed",
                             rid, MAX_ATTEMPTS)
                _record_failure(rid, {"code": "download_lost", "status": 500,
                                      "message": "the worker running this download stopped "
                                                 f"{MAX_ATTEMPTS} times; narrow the request"})
            else:
                _drop_files(d)                     # partial output of the lost run
                write_status(rid, state="queued", done=0, error=None)
                logger.warning("download %s was left running by a stopped worker: queued again",
                               rid)
                queued += 1
        except Exception:                                           # noqa: BLE001
            logger.exception("download %s: cannot recover it", rid)
    return queued


def _strike(rid: str) -> None:
    """A download whose outcome could not be recorded: hold it back, and remove it in the end."""
    n = _strikes[rid] = _strikes.get(rid, 0) + 1
    wait = min(float(settings.DOWNLOADS_POLL_S) * 2 ** n, HOLD_MAX_S)
    _held_until[rid] = time.monotonic() + wait
    if n >= MAX_STRIKES:
        logger.error("download %s failed %d times and its status cannot be written: removed",
                     rid, n)
        remove(rid)
        if result_dir(rid).exists():                  # even that failed: leave it alone for good
            _held_until[rid] = float("inf")
        else:
            _strikes.pop(rid, None)
            _held_until.pop(rid, None)
    else:
        logger.error("download %s: strike %d, next try in %.0f s", rid, n, wait)


def _held() -> set[str]:
    now = time.monotonic()
    return {rid for rid, t in _held_until.items() if t > now}


def run_pending_once(*, stoppable: bool = False) -> int:
    """Take the lock, run every queued download (oldest first). Returns how many were run.

    With ``stoppable`` the call is the runner thread's: it ends when :func:`stop_runner` is called,
    and puts the download it is running back in the queue.
    """
    path = _lock_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    ran = 0
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return 0                                    # another worker is the runner right now
        try:
            recover_orphans()
            while not (stoppable and _stop.is_set()):
                rid = oldest_queued(exclude=_held())
                if rid is None:
                    break
                try:
                    run_download(rid, interruptible=stoppable)
                except Exception:                                   # noqa: BLE001
                    logger.exception("download %s: cannot record its outcome", rid)
                    _strike(rid)
                else:
                    _strikes.pop(rid, None)
                    _held_until.pop(rid, None)
                ran += 1
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
    return ran


def _loop() -> None:
    last_purge = 0.0
    poll = float(settings.DOWNLOADS_POLL_S)
    delay = poll
    while not _stop.is_set():
        try:
            run_pending_once(stoppable=True)
            if time.monotonic() - last_purge > 60:
                purge_expired()
                last_purge = time.monotonic()
            delay = poll
        except Exception:                                           # noqa: BLE001
            logger.exception("downloads runner error")
            delay = min(delay * 2, HOLD_MAX_S)          # a persistent fault (disk) is not retried hot
        _stop.wait(delay)


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
        _thread.join(timeout=STOP_WAIT_S)
    _thread = None
