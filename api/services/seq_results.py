"""On-disk results and job status of /sequences, keyed by ``result_id``.

Layout: ``<root>/<result_id>/{result.fa, manifest.tsv, status.json}``. Files are written
atomically, so any API worker can serve a result or report progress, and an equal request reuses
the finished result. Results older than ``SEQUENCES_RESULT_TTL_S`` are removed on new requests.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
import threading
import time
from collections.abc import Callable
from pathlib import Path

from settings import settings

_ID = re.compile(r"^[0-9a-f]{20}$")
_LIVE: set[str] = set()                    # jobs running in THIS process
_LOCK = threading.Lock()
STALE_AFTER_S = 600


def results_root() -> Path:
    root = Path(settings.SEQUENCES_RESULTS_DIR) if settings.SEQUENCES_RESULTS_DIR else (
        Path(tempfile.gettempdir()) / "annotrieve-sequences")
    root.mkdir(parents=True, exist_ok=True)
    return root


def valid_id(result_id: str) -> bool:
    return bool(_ID.match(result_id or ""))


def result_dir(result_id: str) -> Path:
    if not valid_id(result_id):
        raise ValueError("invalid result id")
    return results_root() / result_id


def write_status(result_id: str, **fields) -> dict:
    d = result_dir(result_id)
    d.mkdir(parents=True, exist_ok=True)
    cur = read_status(result_id, _raw=True) or {}
    cur.update(fields)
    cur["updated"] = time.time()
    tmp = d / "status.json.tmp"
    tmp.write_text(json.dumps(cur))
    os.replace(tmp, d / "status.json")
    return cur


def read_status(result_id: str, *, _raw: bool = False) -> dict | None:
    if not valid_id(result_id):
        return None
    path = result_dir(result_id) / "status.json"
    try:
        status = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    if (not _raw and status.get("state") in ("queued", "running")
            and result_id not in _LIVE and time.time() - status.get("updated", 0) > STALE_AFTER_S):
        status["state"] = "failed"
        status["error"] = {"code": "job_lost", "message": "the worker running this job stopped",
                           "status": 500}
    return status


def purge_expired() -> int:
    cutoff = time.time() - int(settings.SEQUENCES_RESULT_TTL_S)
    removed = 0
    for d in results_root().iterdir():
        if d.is_dir() and valid_id(d.name) and d.name not in _LIVE:
            try:
                if d.stat().st_mtime < cutoff:
                    shutil.rmtree(d, ignore_errors=True)
                    removed += 1
            except OSError:
                pass
    return removed


class Refused(Exception):
    """A request the service will not start now (capacity or per-client limits)."""

    def __init__(self, status: int, code: str, message: str, retry_after: int, **extra) -> None:
        super().__init__(message)
        self.status, self.code, self.retry_after, self.extra = status, code, retry_after, extra


def _entries() -> list[tuple[Path, dict | None]]:
    out = []
    for d in results_root().iterdir():
        if d.is_dir() and valid_id(d.name):
            out.append((d, read_status(d.name)))
    return out


def admit(owner: str) -> None:
    """Check capacity and the client's limits before a new export starts; raises :class:`Refused`.

    Makes room when the store is full by deleting the oldest finished results. Counts come from the
    retained results' ``status.json``, so every worker sees the same numbers. The check is soft: two
    simultaneous requests can both pass it.
    """
    purge_expired()
    now = time.time()
    entries = _entries()

    max_active = int(settings.SEQUENCES_IP_MAX_ACTIVE)
    if max_active > 0:
        active = sum(1 for _d, st in entries if st and st.get("owner") == owner
                     and st.get("state") in ("queued", "running"))
        if active >= max_active:
            raise Refused(429, "too_many_active_exports",
                          f"you already have {active} exports running; wait for one to finish",
                          retry_after=15, active=active, limit=max_active)
    per_hour = int(settings.SEQUENCES_IP_MAX_PER_HOUR)
    if per_hour > 0:
        mine = sorted(st.get("created", 0) for _d, st in entries
                      if st and st.get("owner") == owner and now - st.get("created", 0) < 3600)
        if len(mine) >= per_hour:
            wait = int(max(1, mine[0] + 3600 - now))
            raise Refused(429, "export_rate_limited",
                          f"at most {per_hour} new exports per hour; retry in {wait} s",
                          retry_after=wait, limit_per_hour=per_hour)

    cap = int(settings.SEQUENCES_MAX_RESULTS)
    if cap > 0 and len(entries) >= cap:
        finished = sorted((d for d, st in entries
                           if d.name not in _LIVE and (st is None or st.get("state") in ("done", "failed"))),
                          key=lambda d: d.stat().st_mtime)
        for d in finished[: len(entries) - cap + 1]:
            shutil.rmtree(d, ignore_errors=True)
        if len(entries) - min(len(finished), len(entries) - cap + 1) >= cap:
            raise Refused(503, "sequences_capacity",
                          f"{cap} exports are running; retry shortly", retry_after=30, max_results=cap)


def start_job(result_id: str, work: Callable[[Path, Callable[..., None]], dict], *,
              owner: str = "") -> bool:
    """Run ``work(out_dir, progress)`` in a background thread. False if already running here."""
    with _LOCK:
        if result_id in _LIVE:
            return False
        _LIVE.add(result_id)
    write_status(result_id, state="queued", done=0, error=None, created=time.time(), owner=owner)

    def progress(done, total, counts) -> None:
        write_status(result_id, state="running", done=done, total=total,
                     n_ok=counts.get("ok", 0), n_skipped=counts.get("skipped", 0))

    def run() -> None:
        try:
            write_status(result_id, state="running")
            summary = work(result_dir(result_id), progress)
            write_status(result_id, state="done", **summary)
        except Exception as exc:                                    # noqa: BLE001
            code = getattr(exc, "code", "job_failed")
            write_status(result_id, state="failed", error={
                "code": code, "message": str(exc), "status": getattr(exc, "status", 500)})
        finally:
            with _LOCK:
                _LIVE.discard(result_id)

    threading.Thread(target=run, name=f"sequences-{result_id}", daemon=True).start()
    return True
