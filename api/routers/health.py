"""Health check route: reports what the search API can actually serve."""
from __future__ import annotations

import sqlite3

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from helpers.presence import (
    LOCI_FILE,
    open_request_connections,
    presence_enabled,
    presence_reader,
)

router = APIRouter()


def _corpus_counts(conn: sqlite3.Connection) -> dict | None:
    """How many indexed annotations carry gene names, from the build's ``meta`` table.

    Older builds do not record it; then the key is left out rather than guessed.
    """
    rows = dict(conn.execute(
        "SELECT key, value FROM meta WHERE key IN ('n_annotations', 'n_annotations_with_names')"
    ).fetchall())
    try:
        total = int(rows["n_annotations"])
        named = int(rows["n_annotations_with_names"])
    except (KeyError, TypeError, ValueError):
        return None
    if total <= 0 or named < 0 or named > total:
        return None
    return {"n_annotations": total, "n_annotations_with_names": named}


def _check_presence() -> dict:
    """Open the published index, run a test lookup, and check loci.sqlite matches it."""
    out: dict = {}
    reader = presence_reader()  # 503 HTTPException when nothing is published
    out["data_version"] = reader.data_version
    out["index"] = "ram" if reader.db_path != reader.path else "disk"
    row = reader.conn.execute("SELECT name FROM name_dict LIMIT 1").fetchone()
    out["lookup"] = "ok" if row is not None else "empty"
    corpus = _corpus_counts(reader.conn)
    if corpus is not None:
        out["corpus"] = corpus
    if not (reader.path.parent / LOCI_FILE).is_file():
        out["loci"] = "unpublished"
        return out
    pres, loci, _ = open_request_connections()  # 503 on version mismatch
    try:
        loci.execute("SELECT 1 FROM sqlite_master LIMIT 1").fetchone()
        out["loci"] = "ok"
    finally:
        pres.close()
        loci.close()
    return out


@router.get("/health")
def health():
    if not presence_enabled():
        return {"status": "ok"}
    try:
        checks = _check_presence()
    except Exception as exc:  # any failure means "do not send traffic here"
        detail = getattr(exc, "detail", None)
        reason = (detail.get("code") if isinstance(detail, dict) else None) or type(exc).__name__
        return JSONResponse(status_code=503, content={"status": "unavailable", "reason": reason})
    status = "ok" if checks["lookup"] == "ok" else "degraded"
    return {"status": status, **checks}
