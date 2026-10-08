"""Namespace and accession catalogue routes."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import http_error
from helpers.pagination import DEFAULT_LIMIT, MAX_LIMIT, PAGE_TOKEN_DESC
from services.namespaces import (
    ACCESSION_SORTS,
    get_accession,
    get_namespace,
    list_accessions,
    list_namespaces,
)

router = APIRouter()


@router.get("/namespaces")
def get_namespaces():
    """Namespaces present in the published gene corpus."""
    return list_namespaces(corpus_conn())


@router.get("/namespaces/{namespace:path}/accessions/{accession:path}")
def get_namespace_accession_item(namespace: str, accession: str):
    """One accession summary in a namespace (GO / symbol normalized)."""
    if get_namespace(corpus_conn(), namespace) is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Unknown namespace: {namespace}",
        )
    row = get_accession(corpus_conn(), namespace, accession)
    if row is None:
        raise http_error(
            404,
            "accession_not_found",
            f"Unknown accession: {accession}",
        )
    return row


@router.get("/namespaces/{namespace:path}/accessions")
def get_namespace_accessions(
    namespace: str,
    sort: str = Query("n_genes"),
    sort_order: str = Query("desc"),
    accessions: str | None = Query(
        None,
        description=(
            "Optional comma-separated accession values to keep "
            "(normalized; GO padded to GO:#######)"
        ),
    ),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    next: str | None = Query(None, description=PAGE_TOKEN_DESC),
    previous: str | None = Query(None, description=PAGE_TOKEN_DESC),
):
    """Paginated accessions in one namespace."""
    if sort not in ACCESSION_SORTS:
        raise http_error(
            400,
            "invalid_sort",
            f"sort must be one of: {', '.join(sorted(ACCESSION_SORTS))}",
        )
    if sort_order not in {"asc", "desc"}:
        raise http_error(
            400,
            "invalid_sort_order",
            "sort_order must be asc or desc",
        )
    try:
        page = list_accessions(
            corpus_conn(),
            namespace,
            sort=sort,
            sort_order=sort_order,
            next=next,
            previous=previous,
            limit=limit,
            accessions=accessions,
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
