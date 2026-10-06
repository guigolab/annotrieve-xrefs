"""Opaque pagination cursors for the annotrieve-xrefs API."""
from __future__ import annotations

import base64
import json
from typing import Any, Literal

CursorDir = Literal["next", "prev"]

_CURSOR_VERSION = 1
MAX_CURSOR_CHARS = 4096


class QueryError(ValueError):
    """Client query rejected (bad cursor, sort mismatch, etc.)."""


def encode_cursor(payload: dict[str, Any]) -> str:
    """URL-safe base64 JSON cursor (no padding)."""
    raw = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def decode_cursor(token: str) -> dict[str, Any]:
    """Decode a cursor token; raise QueryError if invalid."""
    if not token or not isinstance(token, str):
        raise QueryError("invalid cursor")
    if len(token) > MAX_CURSOR_CHARS:
        raise QueryError("invalid cursor")
    pad = "=" * (-len(token) % 4)
    try:
        raw = base64.urlsafe_b64decode(token + pad)
        payload = json.loads(raw.decode("utf-8"))
    except (ValueError, json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise QueryError("invalid cursor") from exc
    if not isinstance(payload, dict):
        raise QueryError("invalid cursor")
    if payload.get("v") != _CURSOR_VERSION:
        raise QueryError("invalid cursor")
    return payload


def encode_accessions_cursor(
    *,
    direction: CursorDir,
    sort: str,
    sort_order: str,
    sort_value: int,
    accession: str,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "acc",
            "dir": direction,
            "sort": sort,
            "sort_order": sort_order,
            "n": int(sort_value),
            "a": accession,
        }
    )


def decode_accessions_cursor(
    token: str,
    *,
    sort: str,
    sort_order: str,
) -> tuple[CursorDir, int, str]:
    payload = decode_cursor(token)
    if payload.get("kind") != "acc":
        raise QueryError("invalid cursor")
    if payload.get("sort") != sort or payload.get("sort_order") != sort_order:
        raise QueryError("cursor sort mismatch")
    direction = payload.get("dir")
    sort_value = payload.get("n")
    accession = payload.get("a")
    if (
        direction not in ("next", "prev")
        or not isinstance(sort_value, int)
        or not isinstance(accession, str)
        or len(accession) > 512
    ):
        raise QueryError("invalid cursor")
    return direction, sort_value, accession  # type: ignore[return-value]


def encode_genes_cursor(
    *,
    direction: CursorDir,
    annotation_key: int,
    local_id: int,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "genes",
            "dir": direction,
            "ak": int(annotation_key),
            "lid": int(local_id),
        }
    )


def decode_genes_cursor(token: str) -> tuple[CursorDir, int, int]:
    payload = decode_cursor(token)
    if payload.get("kind") != "genes":
        raise QueryError("invalid cursor")
    direction = payload.get("dir")
    annotation_key = payload.get("ak")
    local_id = payload.get("lid")
    if (
        direction not in ("next", "prev")
        or not isinstance(annotation_key, int)
        or not isinstance(local_id, int)
        or annotation_key < 0
        or local_id < 0
    ):
        raise QueryError("invalid cursor")
    return direction, annotation_key, local_id  # type: ignore[return-value]
