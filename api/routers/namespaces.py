"""Namespace directory, /matches facet, and scoped terms routes."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import http_error
from helpers.pagination import DEFAULT_LIMIT, MAX_LIMIT, MAX_OFFSET, MIN_Q_LEN
from services.namespaces import (
    get_namespace,
    get_term,
    list_namespace_matches,
    list_namespaces,
    list_terms,
)

router = APIRouter()


@router.get("/matches")
def get_matches(
    q: str = Query(..., description="Search string (CURIE, id, symbol, or prose)"),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """
    Namespace facet for *q*: which searchable namespaces match, with counts.

    Every result includes ``match_count``. ``q`` is required.
    """
    try:
        return list_namespace_matches(
            corpus_conn(), q, limit=limit, offset=offset
        )
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=exc.as_detail()) from exc


@router.get("/namespaces")
def get_namespaces(
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """Searchable namespaces present in the published gene corpus."""
    return list_namespaces(corpus_conn(), limit=limit, offset=offset)


@router.get("/namespaces/{namespace:path}/terms/{accession:path}")
def get_namespace_term(namespace: str, accession: str):
    """One term card (CURIE accession) in a searchable namespace."""
    if get_namespace(corpus_conn(), namespace) is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Unknown namespace: {namespace}",
        )
    row = get_term(corpus_conn(), namespace, accession)
    if row is None:
        raise http_error(
            404,
            "term_not_found",
            f"Unknown term: {accession}",
        )
    return row


@router.get("/namespaces/{namespace:path}/terms")
def get_namespace_terms(
    namespace: str,
    q: str | None = Query(
        None,
        description=(
            f"Filter terms (FTS for labeled ns; prefix/exact for identity; "
            f"min {MIN_Q_LEN} chars when searching)"
        ),
    ),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """Paginated terms in one searchable namespace."""
    try:
        page = list_terms(
            corpus_conn(),
            namespace,
            q=q,
            limit=limit,
            offset=offset,
        )
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=exc.as_detail()) from exc
    if page is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Unknown namespace: {namespace}",
        )
    return page


@router.get("/namespaces/{namespace:path}")
def get_namespace_item(namespace: str):
    """One namespace summary."""
    row = get_namespace(corpus_conn(), namespace)
    if row is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Unknown namespace: {namespace}",
        )
    return row
