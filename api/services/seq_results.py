"""On-disk downloads (``/downloads``): request, status, queue order, daily quota and clean-up.

Layout: ``<root>/<id>/{request.json, status.json, result.fa.gz, manifest.tsv}`` and
``<root>/.quota/<owner>.json``. ``id`` identifies a request (normalised body + data version), so an
equal request reuses one download. Files are written atomically, so every API worker sees the same
state. The queue is the set of ``queued`` statuses ordered by ``seq``; one runner at a time takes
the oldest (see ``seq_runner``).

How long files stay: ``DOWNLOAD_KEEP_AFTER_FETCH_S`` after the first complete download,
``DOWNLOAD_KEEP_UNFETCHED_S`` if never fetched, ``DOWNLOAD_KEEP_FAILED_S`` for a failed one.
"""
from __future__ import annotations

import fcntl
import json
import os
import re
import shutil
import tempfile
import threading
import time
from collections.abc import Collection, Iterator
from contextlib import contextmanager
from pathlib import Path

from settings import settings

_ID = re.compile(r"^[0-9a-f]{20}$")
DAY_S = 24 * 3600


def stale_after_s() -> float:
    """A running download that has not updated its status for this long is lost."""
    return float(settings.SEQUENCES_DEADLINE_S) + 60.0


def results_root() -> Path:
    root = Path(settings.SEQUENCES_RESULTS_DIR) if settings.SEQUENCES_RESULTS_DIR else (
        Path(tempfile.gettempdir()) / "annotrieve-sequences")
    root.mkdir(parents=True, exist_ok=True)
    return root


def valid_id(result_id: str) -> bool:
    return bool(_ID.match(result_id or ""))


def result_dir(result_id: str) -> Path:
    if not valid_id(result_id):
        raise ValueError("invalid download id")
    return results_root() / result_id


def _write_atomic(path: Path, text: str) -> None:
    """Write ``path`` through a temp file of its own (two writers never share one) and rename it."""
    tmp = path.with_name(f"{path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    try:
        tmp.write_text(text)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def write_status(result_id: str, **fields) -> dict:
    d = result_dir(result_id)
    d.mkdir(parents=True, exist_ok=True)
    cur = read_status(result_id, _raw=True) or {}
    cur.update(fields)
    cur["updated"] = time.time()
    _write_atomic(d / "status.json", json.dumps(cur))
    return cur


def read_status(result_id: str, *, _raw: bool = False) -> dict | None:
    if not valid_id(result_id):
        return None
    try:
        status = json.loads((result_dir(result_id) / "status.json").read_text())
    except (OSError, ValueError):
        return None
    if (not _raw and status.get("state") == "running"
            and time.time() - status.get("updated", 0) > stale_after_s()):
        status["state"] = "failed"
        status["error"] = {"code": "download_lost", "message": "the worker running this download "
                           "stopped", "status": 500}
    return status


def fasta_path(result_id: str) -> Path:
    return result_dir(result_id) / "result.fa.gz"


def remove(result_id: str) -> None:
    shutil.rmtree(result_dir(result_id), ignore_errors=True)


def mark_downloaded(result_id: str) -> None:
    """Called once a client has received the whole file: starts the short keep-time."""
    if read_status(result_id, _raw=True) is not None:
        write_status(result_id, downloaded_at=time.time())


# --- daily quota (rolling 24 h per client; one small file per owner) --------------------------------

def _quota_file(owner: str) -> Path:
    d = results_root() / ".quota"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{re.sub(r'[^0-9a-f]', '', owner)[:32] or 'unknown'}.json"


def _quota_times(owner: str, now: float) -> list[float]:
    try:
        raw = json.loads(_quota_file(owner).read_text())
    except (OSError, ValueError):
        return []
    return sorted(t for t in raw if isinstance(t, (int, float)) and now - t < DAY_S)


def quota_status(owner: str) -> dict:
    """Used and remaining downloads in the last 24 h, and when the oldest one stops counting."""
    now = time.time()
    times = _quota_times(owner, now)
    limit = int(settings.DOWNLOADS_PER_DAY)
    return {"limit": limit, "used": len(times),
            "remaining": max(limit - len(times), 0) if limit > 0 else None,
            "window_hours": 24,
            "resets_at": (times[0] + DAY_S) if times and limit > 0 and len(times) >= limit else None}


def _quota_take(owner: str) -> None:
    now = time.time()
    times = _quota_times(owner, now) + [now]
    _write_atomic(_quota_file(owner), json.dumps(times))


# --- clean-up, capacity and admission ----------------------------------------------------------------

def purge_expired() -> int:
    now = time.time()
    removed = 0
    for d in results_root().iterdir():
        if not (d.is_dir() and valid_id(d.name)):
            continue
        st = read_status(d.name, _raw=True)
        if st is None:                       # directory without status: leftover of a crash
            try:
                if now - d.stat().st_mtime > 3600:
                    shutil.rmtree(d, ignore_errors=True)
                    removed += 1
            except OSError:
                pass
            continue
        state, upd = st.get("state"), st.get("updated", 0)
        if state == "done":
            if st.get("downloaded_at"):
                old = now - st["downloaded_at"] > int(settings.DOWNLOAD_KEEP_AFTER_FETCH_S)
            else:
                old = now - upd > int(settings.DOWNLOAD_KEEP_UNFETCHED_S)
        elif state == "failed":
            old = now - upd > int(settings.DOWNLOAD_KEEP_FAILED_S)
        else:                                # queued or running: never expire (lost ones show failed)
            old = state == "running" and now - upd > stale_after_s() + int(settings.DOWNLOAD_KEEP_FAILED_S)
        if old:
            shutil.rmtree(d, ignore_errors=True)
            removed += 1
    return removed


class Refused(Exception):
    """A request the service will not start now (capacity or per-client limits)."""

    def __init__(self, status: int, code: str, message: str, retry_after: int, **extra) -> None:
        super().__init__(message)
        self.status, self.code, self.retry_after, self.extra = status, code, retry_after, extra


def entries() -> list[tuple[Path, dict | None]]:
    out = []
    for d in results_root().iterdir():
        if d.is_dir() and valid_id(d.name):
            out.append((d, read_status(d.name)))
    return out


def queue_position(result_id: str) -> int | None:
    """1 = next to run (or running). None when the download is not waiting."""
    me = read_status(result_id)
    if not me or me.get("state") != "queued":
        return None
    ahead = sum(1 for _d, st in entries() if st and st.get("state") == "queued"
                and st.get("seq", 0) < me.get("seq", 0))
    return ahead + 1


@contextmanager
def admission_lock() -> Iterator[None]:
    """One admission at a time for the whole container (every worker).

    The per-client limits and the queue size are counted from the files on disk, so two requests
    that are checked at the same moment would both pass. The lock file is on the container's own
    disk, next to the runner's lock (``flock`` is unreliable on the network drive).
    """
    raw = (settings.DOWNLOADS_LOCK_FILE or "").strip()
    base = Path(raw) if raw else Path(tempfile.gettempdir()) / "annotrieve-downloads.lock"
    path = base.with_name(base.name + ".admit")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def admit(owner: str) -> None:
    """Check the client's limits and the capacity before a new download; raises :class:`Refused`.

    A download that is admitted uses one place of the client's daily quota. Call it under
    :func:`admission_lock`, together with :func:`enqueue`, so the counts it reads stay true.
    """
    purge_expired()
    ents = entries()

    max_active = int(settings.SEQUENCES_IP_MAX_ACTIVE)
    if max_active > 0:
        active = sum(1 for _d, st in ents if st and st.get("owner") == owner
                     and st.get("state") in ("queued", "running"))
        if active >= max_active:
            raise Refused(429, "too_many_active_downloads",
                          f"you already have {active} downloads waiting or running; wait for one "
                          "to finish", retry_after=30, active=active, limit=max_active)
    q = quota_status(owner)
    if q["limit"] > 0 and q["used"] >= q["limit"]:
        wait = int(max(1, (q["resets_at"] or time.time() + DAY_S) - time.time()))
        raise Refused(429, "daily_quota_exceeded",
                      f"at most {q['limit']} new downloads per 24 hours per client; the next one is "
                      f"free in {wait} s", retry_after=wait, quota=q)

    max_queued = int(settings.DOWNLOADS_MAX_QUEUED)
    if max_queued > 0:
        queued = sum(1 for _d, st in ents if st and st.get("state") == "queued")
        if queued >= max_queued:
            raise Refused(503, "downloads_queue_full",
                          f"{queued} downloads are waiting; retry later", retry_after=300,
                          queued=queued, max_queued=max_queued)
    cap = int(settings.SEQUENCES_MAX_RESULTS)
    if cap > 0 and len(ents) >= cap:
        raise Refused(503, "downloads_capacity",
                      f"{cap} downloads are stored; retry once some have been fetched or expired",
                      retry_after=300, max_results=cap)
    _quota_take(owner)


def enqueue(result_id: str, body: dict, *, owner: str, total: int, data_version: str) -> dict:
    """Store the request and put the download at the end of the queue."""
    d = result_dir(result_id)
    d.mkdir(parents=True, exist_ok=True)
    _write_atomic(d / "request.json", json.dumps(body))
    return write_status(result_id, state="queued", seq=time.time_ns(), created=time.time(),
                        owner=owner, done=0, total=total, data_version=data_version, error=None,
                        downloaded_at=None, attempts=0)


def running_ids() -> list[str]:
    """Downloads whose status says ``running``, as written (not the ``download_lost`` view).

    Only the runner holding the downloads lock may call this to find orphans: while it holds the
    lock, no other process runs a download, so a ``running`` status is the leftover of a worker
    that stopped.
    """
    out = []
    for d in results_root().iterdir():
        if d.is_dir() and valid_id(d.name):
            st = read_status(d.name, _raw=True)
            if st and st.get("state") == "running":
                out.append(d.name)
    return out


def oldest_queued(exclude: Collection[str] = ()) -> str | None:
    """The queued download with the lowest ``seq``, skipping the ids in ``exclude``."""
    best: tuple[int, str] | None = None
    for d, st in entries():
        if st and st.get("state") == "queued" and d.name not in exclude:
            key = (int(st.get("seq", 0)), d.name)
            if best is None or key < best:
                best = key
    return best[1] if best else None
