"""Per-annotation gene and xref routes (shard-backed)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Path, Query

from helpers.cursor import QueryError
from helpers.deps import annotation_shard
from helpers.errors import error_detail, http_error
from helpers.match import MATCH_ANY
from helpers.pagination import (
    DEFAULT_LIMIT,
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

_MATCH_DESC = "any (union over terms, default) or all (intersection over terms)"

# Validation-style QueryError codes → HTTP 422 (align with /hits)
_VALIDATION_CODES = frozenset(
    {
        "names_required",
        "too_many_names",
        "prefix_too_short",
        "prefix_too_broad",
        "invalid_match",
    }
)


def _http_for_query_error(exc: QueryError) -> HTTPException:
    status = 422 if exc.code in _VALIDATION_CODES else 400
    return HTTPException(status_code=status, detail=exc.as_detail())


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
    names: str | None = Query(
        None,
        description=(
            "Comma-separated gene names; optional trailing * for prefix "
            "(symbol∪alias)"
        ),
    ),
    match: str = Query(MATCH_ANY, description=_MATCH_DESC),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
):
    """Paginated gene cards for one annotation (optional names filter)."""
    conn = annotation_shard(annotation_id)
    try:
        try:
            return list_annotation_genes(
                conn,
                annotation_id,
                names=names,
                match=match,
                limit=limit,
                offset=offset,
            )
        except QueryError as exc:
            raise _http_for_query_error(exc) from exc
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
        row = get_gene(conn, local_id)
    finally:
        conn.close()
    if row is None:
        raise http_error(
            404,
            "gene_not_found",
            f"Gene local_id {local_id} not found in {annotation_id}",
        )
    return row
