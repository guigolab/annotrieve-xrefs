"""Per-annotation gene and xref routes (shard-backed)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Path, Query

from helpers.cursor import QueryError
from helpers.deps import annotation_shard
from helpers.errors import http_error
from helpers.pagination import DEFAULT_LIMIT, MAX_LIMIT, PAGE_TOKEN_DESC
from services.annotations import get_gene, list_annotation_genes, list_gene_xrefs

router = APIRouter()


@router.get("/annotations/{annotation_id}/genes")
def get_annotation_genes(
    annotation_id: str,
    q: str | None = Query(
        None,
        description=(
            "Optional search: a CURIE (contains ':') filters by xref; "
            "otherwise casefold prefix match on primary_name"
        ),
    ),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    next: str | None = Query(None, description=PAGE_TOKEN_DESC),
    previous: str | None = Query(None, description=PAGE_TOKEN_DESC),
):
    """Paginated gene cards for one annotation (optional q filter)."""
    conn = annotation_shard(annotation_id)
    try:
        try:
            return list_annotation_genes(
                conn,
                annotation_id,
                q=q,
                next=next,
                previous=previous,
                limit=limit,
            )
        except QueryError as exc:
            raise HTTPException(status_code=400, detail=exc.as_detail()) from exc
    finally:
        conn.close()


@router.get("/annotations/{annotation_id}/genes/{local_id}/xrefs")
def get_annotation_gene_xrefs(
    annotation_id: str,
    local_id: int = Path(..., ge=0),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    next: str | None = Query(None, description=PAGE_TOKEN_DESC),
    previous: str | None = Query(None, description=PAGE_TOKEN_DESC),
):
    """Paginated xref chips for one gene in this annotation's shard."""
    conn = annotation_shard(annotation_id)
    try:
        try:
            page = list_gene_xrefs(
                conn,
                annotation_id,
                local_id,
                next=next,
                previous=previous,
                limit=limit,
            )
        except QueryError as exc:
            raise HTTPException(status_code=400, detail=exc.as_detail()) from exc
    finally:
        conn.close()
    if page is None:
        raise http_error(
            404,
            "gene_not_found",
            f"Gene local_id {local_id} not found in {annotation_id}",
        )
    return page


@router.get("/annotations/{annotation_id}/genes/{local_id}")
def get_annotation_gene(
    annotation_id: str,
    local_id: int = Path(..., ge=0),
):
    """One gene locus card from this annotation's gene corpus shard."""
    conn = annotation_shard(annotation_id)
    try:
        row = get_gene(conn, local_id, annotation_id=annotation_id)
    finally:
        conn.close()
    if row is None:
        raise http_error(
            404,
            "gene_not_found",
            f"Gene local_id {local_id} not found in {annotation_id}",
        )
    return row
