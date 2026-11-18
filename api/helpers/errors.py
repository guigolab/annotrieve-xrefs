"""Shared HTTP error body shape for the annotrieve-xrefs API."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # fastapi is imported lazily so offline tools can use the services
    from fastapi import HTTPException


def error_detail(code: str, message: str, **extra: Any) -> dict[str, Any]:
    """Structured error payload: ``code``, ``message``, plus optional extras."""
    detail: dict[str, Any] = {"code": code, "message": message}
    detail.update(extra)
    return detail


def http_error(
    status_code: int,
    code: str,
    message: str,
    **extra: Any,
) -> HTTPException:
    """``HTTPException`` whose ``detail`` is an ``error_detail`` object."""
    from fastapi import HTTPException

    return HTTPException(
        status_code=status_code,
        detail=error_detail(code, message, **extra),
    )
