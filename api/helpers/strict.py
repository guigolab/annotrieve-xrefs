"""Reject query parameters an endpoint does not know, instead of silently ignoring them.

FastAPI drops unknown query parameters. A client that asks for ``format=tsv`` or ``min_matched=2``
on an endpoint that lacks them then gets an unfiltered JSON answer that looks right. This turns that
into a 422 listing the parameters the endpoint accepts.
"""
from __future__ import annotations

from collections.abc import Iterable

from fastapi import Request

from helpers.errors import http_error


def reject_unknown_params(request: Request, allowed: Iterable[str]) -> None:
    allowed_set = set(allowed)
    unknown = sorted(set(request.query_params) - allowed_set)
    if unknown:
        raise http_error(
            422, "unknown_parameter",
            f"unknown query parameter(s): {', '.join(unknown)}",
            unknown=unknown, allowed=sorted(allowed_set),
        )
