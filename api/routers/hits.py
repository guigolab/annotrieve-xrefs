"""Cross-annotation hits and annotation-collapsed search routes."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from helpers.concurrency import BusyError, hits_slot
from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import error_detail
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_CURIES_GET,
    MAX_CURIES_POST,
    MAX_LIMIT,
    PAGE_TOKEN_DESC,
)
from services.hits import MATCH_ANY, list_hit_annotations, list_hits

router = APIRouter()

_MATCH_DESC = "any (union, default) or all (intersection)"


class HitsSearchBody(BaseModel):
    """JSON body for ``POST /hits`` and ``POST /hits/annotations``."""

    curies: list[str]
    taxid: int | None = Field(default=None, ge=1)
    match: str = Field(default=MATCH_ANY, description=_MATCH_DESC)
    limit: int = Field(default=DEFAULT_LIMIT, ge=1, le=MAX_LIMIT)
    next: str | None = None
    previous: str | None = None


def _busy_response() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={
            "detail": error_detail(
                "hits_busy",
                "too many concurrent hits requests; retry shortly",
            )
        },
        headers={"Retry-After": "1"},
    )


def _hits_response(
    curies: str | list[str],
    *,
    taxid: int | None,
    match: str,
    next: str | None,
    previous: str | None,
    limit: int,
    max_curies: int,
):
    try:
        with hits_slot():
            return list_hits(
                corpus_conn(),
                curies,
                taxid=taxid,
                match=match,
                next=next,
                previous=previous,
                limit=limit,
                max_curies=max_curies,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=exc.as_detail()) from exc


def _hit_annotations_response(
    curies: str | list[str],
    *,
    taxid: int | None,
    match: str,
    next: str | None,
    previous: str | None,
    limit: int,
    max_curies: int,
):
    try:
        with hits_slot():
            return list_hit_annotations(
                corpus_conn(),
                curies,
                taxid=taxid,
                match=match,
                next=next,
                previous=previous,
                limit=limit,
                max_curies=max_curies,
            )
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=exc.as_detail()) from exc


@router.get("/hits/annotations")
def get_hit_annotations(
    curies: str = Query(
        ...,
        description="Comma-separated CURIEs, e.g. symbol:tp53,GO:0008150",
    ),
    taxid: int | None = Query(
        None,
        ge=1,
        description="Optional NCBI taxid (species or ancestor via annotation_lineage)",
    ),
    match: str = Query(MATCH_ANY, description=_MATCH_DESC),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    next: str | None = Query(None, description=PAGE_TOKEN_DESC),
    previous: str | None = Query(None, description=PAGE_TOKEN_DESC),
):
    """Annotations matching the given CURIEs (one row per annotation)."""
    return _hit_annotations_response(
        curies,
        taxid=taxid,
        match=match,
        next=next,
        previous=previous,
        limit=limit,
        max_curies=MAX_CURIES_GET,
    )


@router.post("/hits/annotations")
def post_hit_annotations(body: HitsSearchBody):
    """Annotation-collapsed hits for a larger CURIE list (JSON body)."""
    return _hit_annotations_response(
        body.curies,
        taxid=body.taxid,
        match=body.match,
        next=body.next,
        previous=body.previous,
        limit=body.limit,
        max_curies=MAX_CURIES_POST,
    )


@router.get("/hits")
def get_hits(
    curies: str = Query(
        ...,
        description="Comma-separated CURIEs, e.g. symbol:tp53,GO:0008150",
    ),
    taxid: int | None = Query(
        None,
        ge=1,
        description="Optional NCBI taxid (species or ancestor via annotation_lineage)",
    ),
    match: str = Query(MATCH_ANY, description=_MATCH_DESC),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    next: str | None = Query(None, description=PAGE_TOKEN_DESC),
    previous: str | None = Query(None, description=PAGE_TOKEN_DESC),
):
    """Cross-annotation hits matching the given CURIEs (keyset pagination)."""
    return _hits_response(
        curies,
        taxid=taxid,
        match=match,
        next=next,
        previous=previous,
        limit=limit,
        max_curies=MAX_CURIES_GET,
    )


@router.post("/hits")
def post_hits(body: HitsSearchBody):
    """Hits for a larger CURIE list (JSON body, keyset pagination)."""
    return _hits_response(
        body.curies,
        taxid=body.taxid,
        match=body.match,
        next=body.next,
        previous=body.previous,
        limit=body.limit,
        max_curies=MAX_CURIES_POST,
    )
