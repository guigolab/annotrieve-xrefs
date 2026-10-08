"""Opaque pagination cursors for the annotrieve-xrefs API."""
from __future__ import annotations

import base64
import hashlib
import json
from typing import Any, Literal

from helpers.errors import error_detail

CursorDir = Literal["next", "prev"]

_CURSOR_VERSION = 1
MAX_CURSOR_CHARS = 4096


_MESSAGE_CODES = {
    "invalid cursor": "invalid_cursor",
    "cursor sort mismatch": "cursor_sort_mismatch",
    "page token does not match current query filters": "cursor_filter_mismatch",
}


class QueryError(ValueError):
    """Client query rejected (bad cursor, sort mismatch, etc.)."""

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        **extra: Any,
    ) -> None:
        super().__init__(message)
        self.code = code if code is not None else _MESSAGE_CODES.get(
            message, "invalid_query"
        )
        self.extra = extra

    def as_detail(self) -> dict[str, Any]:
        """Structured 4xx body: ``code``, ``message``, plus any extras."""
        return error_detail(self.code, str(self), **self.extra)


def filter_fingerprint(parts: dict[str, Any]) -> str:
    """Stable 16-char fingerprint of canonical filter JSON."""
    raw = json.dumps(parts, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:16]


def assert_cursor_filter(payload: dict[str, Any], expected_f: str) -> None:
    """Raise cursor_filter_mismatch when the token's ``f`` does not match."""
    got = payload.get("f")
    if got != expected_f:
        raise QueryError(
            "page token does not match current query filters",
            code="cursor_filter_mismatch",
        )


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


def resolve_page_token(
    *,
    next: str | None = None,
    previous: str | None = None,
) -> tuple[CursorDir | None, str | None]:
    """
    Map mutually exclusive ``next`` / ``previous`` query values to a page seek.

    Returns ``(direction, token)`` or ``(None, None)`` when neither is set.
    """
    has_next = bool(next)
    has_prev = bool(previous)
    if has_next and has_prev:
        raise QueryError(
            "pass only one of next or previous",
            code="conflicting_pagination",
        )
    if has_next:
        return "next", next
    if has_prev:
        return "prev", previous
    return None, None


def encode_accessions_cursor(
    *,
    sort: str,
    sort_order: str,
    sort_value: int,
    accession: str,
    filter_f: str,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "acc",
            "sort": sort,
            "sort_order": sort_order,
            "n": int(sort_value),
            "a": accession,
            "f": filter_f,
        }
    )


def decode_accessions_cursor(
    token: str,
    *,
    sort: str,
    sort_order: str,
    filter_f: str,
) -> tuple[int, str]:
    payload = decode_cursor(token)
    if payload.get("kind") != "acc":
        raise QueryError("invalid cursor")
    if payload.get("sort") != sort or payload.get("sort_order") != sort_order:
        raise QueryError(
            "page token does not match current query filters",
            code="cursor_filter_mismatch",
        )
    assert_cursor_filter(payload, filter_f)
    sort_value = payload.get("n")
    accession = payload.get("a")
    if (
        not isinstance(sort_value, int)
        or isinstance(sort_value, bool)
        or not isinstance(accession, str)
        or len(accession) > 512
    ):
        raise QueryError("invalid cursor")
    return sort_value, accession


def encode_hits_cursor(
    *,
    annotation_key: int,
    local_id: int,
    filter_f: str,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "hits",
            "ak": int(annotation_key),
            "lid": int(local_id),
            "f": filter_f,
        }
    )


def decode_hits_cursor(token: str, *, filter_f: str) -> tuple[int, int]:
    payload = decode_cursor(token)
    if payload.get("kind") != "hits":
        raise QueryError("invalid cursor")
    assert_cursor_filter(payload, filter_f)
    annotation_key = payload.get("ak")
    local_id = payload.get("lid")
    if (
        not isinstance(annotation_key, int)
        or isinstance(annotation_key, bool)
        or not isinstance(local_id, int)
        or isinstance(local_id, bool)
        or annotation_key < 0
        or local_id < 0
    ):
        raise QueryError("invalid cursor")
    return annotation_key, local_id


def encode_hits_annotations_cursor(
    *,
    annotation_key: int,
    filter_f: str,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "hits_ann",
            "ak": int(annotation_key),
            "f": filter_f,
        }
    )


def decode_hits_annotations_cursor(token: str, *, filter_f: str) -> int:
    payload = decode_cursor(token)
    if payload.get("kind") != "hits_ann":
        raise QueryError("invalid cursor")
    assert_cursor_filter(payload, filter_f)
    annotation_key = payload.get("ak")
    if (
        not isinstance(annotation_key, int)
        or isinstance(annotation_key, bool)
        or annotation_key < 0
    ):
        raise QueryError("invalid cursor")
    return annotation_key


def encode_ann_genes_cursor(*, local_id: int, filter_f: str) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "ann_genes",
            "lid": int(local_id),
            "f": filter_f,
        }
    )


def decode_ann_genes_cursor(token: str, *, filter_f: str) -> int:
    payload = decode_cursor(token)
    if payload.get("kind") != "ann_genes":
        raise QueryError("invalid cursor")
    assert_cursor_filter(payload, filter_f)
    local_id = payload.get("lid")
    if (
        not isinstance(local_id, int)
        or isinstance(local_id, bool)
        or local_id < 0
    ):
        raise QueryError("invalid cursor")
    return local_id


def encode_gene_xrefs_cursor(
    *,
    namespace: str,
    accession: str,
) -> str:
    return encode_cursor(
        {
            "v": _CURSOR_VERSION,
            "kind": "gene_xrefs",
            "ns": namespace,
            "a": accession,
        }
    )


def decode_gene_xrefs_cursor(token: str) -> tuple[str, str]:
    payload = decode_cursor(token)
    if payload.get("kind") != "gene_xrefs":
        raise QueryError("invalid cursor")
    namespace = payload.get("ns")
    accession = payload.get("a")
    if (
        not isinstance(namespace, str)
        or not isinstance(accession, str)
        or not namespace
        or not accession
        or len(namespace) > 512
        or len(accession) > 512
    ):
        raise QueryError("invalid cursor")
    return namespace, accession
