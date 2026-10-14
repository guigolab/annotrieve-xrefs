"""Per-annotation gene and xref routes (shard-backed)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Path, Query

from helpers.cursor import QueryError
from helpers.deps import annotation_shard
from helpers.errors import error_detail, http_error
from helpers.match import MATCH_ANY
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_CURIES_GET,
    MAX_LIMIT,
    MAX_OFFSET,
)
from services.annotations import (
    get_gene,
    list_annotation_genes,
    list_annotation_namespaces,
    list_gene_xrefs,
)

router = APIRouter()


@router.get("/annotations/{annotation_id}/namespaces")
def get_annotation_namespaces(annotation_id: str):
    """Namespaces present in this annotation's shard (all gene_xref ns)."""
    conn = annotation_shard(annotation_id)
    try:
        return list_annotation_namespaces(conn, annotation_id)
    finally:
        conn.close()


@router.get("/annotations/{annotation_id}/genes")
def get_annotation_genes(
    annotation_id: str,
    q: str | None = Query(
        None,
        description=(
            "Optional casefold prefix match on primary_name (max 64 chars)"
        ),
    ),
    curies: str | None = Query(
        None,
        description="Comma-separated CURIEs, e.g. symbol:tp53,GO:0008150",
    ),
    match: str = Query(
        MATCH_ANY,
        description="any (union, default) or all (intersection); only with curies",
    ),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """Paginated gene cards for one annotation (optional q / curies filters)."""
    conn = annotation_shard(annotation_id)
    try:
        try:
            return list_annotation_genes(
                conn,
                annotation_id,
                q=q,
                curies=curies,
                match=match,
                limit=limit,
                offset=offset,
                max_curies=MAX_CURIES_GET,
            )
        except QueryError as exc:
            raise HTTPException(status_code=400, detail=exc.as_detail()) from exc
        except ValueError as exc:
            raise HTTPException(
                status_code=400,
                detail=error_detail("invalid_offset", str(exc)),
            ) from exc
    finally:
        conn.close()


@router.get("/annotations/{annotation_id}/genes/{local_id}/xrefs")
def get_annotation_gene_xrefs(
    annotation_id: str,
    local_id: int = Path(..., ge=0),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """Paginated xref chips for one gene in this annotation's shard."""
    conn = annotation_shard(annotation_id)
    try:
        try:
            page = list_gene_xrefs(
                conn,
                annotation_id,
                local_id,
                limit=limit,
                offset=offset,
            )
        except ValueError as exc:
            raise HTTPException(
                status_code=400,
                detail=error_detail("invalid_offset", str(exc)),
            ) from exc
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
