"""Query rejection helpers for the annotrieve-xrefs API."""
from __future__ import annotations

from typing import Any

from helpers.errors import error_detail


_MESSAGE_CODES = {
    "invalid cursor": "invalid_cursor",
    "cursor sort mismatch": "cursor_sort_mismatch",
    "page token does not match current query filters": "cursor_filter_mismatch",
}


class QueryError(ValueError):
    """Client query rejected (bad filters, pagination, etc.)."""

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
