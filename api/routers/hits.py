"""Gene-name presence search (``/hits``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from helpers.concurrency import BusyError, hits_slot
from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import error_detail
from helpers.pagination import MAX_NAMES_GET, MAX_NAMES_POST
from services.hits import list_name_hits

router = APIRouter()

# Validation-style QueryError codes → HTTP 422
_VALIDATION_CODES = frozenset(
    {
        "names_required",
        "too_many_names",
        "too_many_taxids",
        "wildcard_not_supported",
        "invalid_taxid",
    }
)


class HitsSearchBody(BaseModel):
    """JSON body for ``POST /hits``.

    Cap is enforced in ``parse_name_list`` (casefold-deduped) so GET and POST
    share ``too_many_names`` — do not use list ``max_length`` here.
    """

    names: list[str] = Field(...)
    include_taxids: list[int] | None = None
    exclude_taxids: list[int] | None = None


def _busy_response() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={
            "detail": error_detail(
                "hits_busy",
                "too many concurrent search requests; retry shortly",
            )
        },
        headers={"Retry-After": "1"},
    )


def _http_for_query_error(exc: QueryError) -> HTTPException:
    status = 422 if exc.code in _VALIDATION_CODES else 400
    return HTTPException(status_code=status, detail=exc.as_detail())


@router.get("/hits")
def get_hits(
    names: str = Query(
        ...,
        description=(
            "Comma-separated exact gene names (symbol∪alias). "
            "Wildcards are not supported; use GET /names?prefix=."
        ),
    ),
    include_taxids: str | None = Query(
        None,
        description="Comma-separated NCBI taxids to include (lineage OR)",
    ),
    exclude_taxids: str | None = Query(
        None,
        description="Comma-separated NCBI taxids to exclude (lineage)",
    ),
):
    """
    Presence search by exact gene symbol/synonym names.

    Compact envelope: ``terms``, ``fields``/``rows`` with ``matched``
    (term indices), ``release_date``, and ``assembly_name``.
    """
    try:
        with hits_slot():
            return list_name_hits(
                corpus_conn(),
                names,
                include_taxids=include_taxids,
                exclude_taxids=exclude_taxids,
                max_names=MAX_NAMES_GET,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc


@router.post("/hits")
def post_hits(body: HitsSearchBody):
    """Presence search for a larger name list (JSON body). Prefer GET when ≤200."""
    try:
        with hits_slot():
            return list_name_hits(
                corpus_conn(),
                body.names,
                include_taxids=body.include_taxids,
                exclude_taxids=body.exclude_taxids,
                max_names=MAX_NAMES_POST,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc
