"""Gene-name presence search (``/hits``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from helpers.concurrency import BusyError, hits_slot
from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import error_detail
from helpers.match import MATCH_ANY
from helpers.pagination import MAX_NAMES_GET, MAX_NAMES_POST
from services.hits import list_name_hits

router = APIRouter()

_MATCH_DESC = "any (union over terms, default) or all (intersection over terms)"

# Validation-style QueryError codes → HTTP 422
_VALIDATION_CODES = frozenset(
    {
        "names_required",
        "too_many_names",
        "too_many_taxids",
        "prefix_too_short",
        "prefix_too_broad",
        "invalid_match",
        "invalid_taxid",
    }
)


class HitsSearchBody(BaseModel):
    """JSON body for ``POST /hits``."""

    names: list[str]
    include_taxids: list[int] | None = None
    exclude_taxids: list[int] | None = None
    match: str = Field(default=MATCH_ANY, description=_MATCH_DESC)


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
            "Comma-separated gene names; optional trailing * for prefix "
            "(symbol∪alias)"
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
    match: str = Query(MATCH_ANY, description=_MATCH_DESC),
):
    """
    Presence search by gene symbol/synonym names.

    Compact envelope: ``terms``, ``keys``, ``fields``/``rows`` with
    ``n_matched`` / ``matched`` (indices into ``keys``),
    ``release_date``, and ``assembly_name``.
    """
    try:
        with hits_slot():
            return list_name_hits(
                corpus_conn(),
                names,
                include_taxids=include_taxids,
                exclude_taxids=exclude_taxids,
                match=match,
                max_names=MAX_NAMES_GET,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc


@router.post("/hits")
def post_hits(body: HitsSearchBody):
    """Presence search for a larger name list (JSON body)."""
    try:
        with hits_slot():
            return list_name_hits(
                corpus_conn(),
                body.names,
                include_taxids=body.include_taxids,
                exclude_taxids=body.exclude_taxids,
                match=body.match,
                max_names=MAX_NAMES_POST,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc
