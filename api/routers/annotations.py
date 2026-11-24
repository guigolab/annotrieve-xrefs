"""Per-annotation gene and xref routes (shard-backed)."""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Path, Query

from helpers.concurrency import annotations_slot
from helpers.cursor import QueryError
from helpers.deps import annotation_shard
from helpers.errors import error_detail, http_error, http_for_query_error
from helpers.match import MATCH_ANY
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_LIMIT,
    MAX_OFFSET,
    resolve_offset,
)

from services.annotations import (
    get_gene,
    list_annotation_genes,
    list_annotation_namespaces,
    list_gene_xrefs,
    list_namespace_accessions,
    resolve_path_namespace,
)

router = APIRouter(dependencies=[Depends(annotations_slot)])

_MATCH_DESC = "any (union over terms, default) or all (intersection over terms)"



def _one_prefix(q: str | None, prefix: str | None) -> str | None:
    """``q`` and ``prefix`` are the same filter; giving both with different values is an error."""
    if q and prefix and q != prefix:
        raise QueryError("q and prefix are the same filter: give one", code="invalid_combination",
                         q=q, prefix=prefix)
    return prefix or q


@router.get("/annotations/{annotation_id}/namespaces")
def get_annotation_namespaces(annotation_id: str):
    """Namespaces present in this annotation's shard, with the number of distinct accessions.

    ``results`` lists everything at once (it is never paged): ``limit`` and ``offset`` in the answer
    only repeat ``total`` and 0, and ``next_cursor`` is always null. Namespaces are the places an
    accession can come from; use one in the accessions route below to browse it. ``symbol`` and
    ``alias`` are not databases: they are the two kinds of gene name that ``/hits`` searches.
    """
    conn = annotation_shard(annotation_id)
    try:
        return list_annotation_namespaces(conn, annotation_id)
    finally:
        conn.close()


@router.get(
    "/annotations/{annotation_id}/namespaces/{namespace:path}/accessions"
)
def get_annotation_namespace_accessions(
    annotation_id: str,
    namespace: str = Path(
        ...,
        description=(
            "Stored namespace or public short prefix (swiss-prot). "
            "Slash forms: literal path or percent-encoded (%2F)."
        ),
    ),
    q: str | None = Query(
        None,
        description=(
            "Optional accession prefix (min 3 characters when set); "
            "omit to browse all. Alias of `prefix`."
        ),
    ),
    prefix: str | None = Query(
        None, description="Optional accession prefix (min 3 characters); same as `q`"),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
    cursor: str | None = Query(
        None, description="`next_cursor` of the previous page; use it instead of `offset`"),
):
    """Paginated accessions under one namespace on this annotation's shard."""
    try:
        q = _one_prefix(q, prefix)
        offset = resolve_offset(cursor, offset)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
    stored_ns = resolve_path_namespace(namespace)
    if stored_ns is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Namespace {namespace!r} not found in {annotation_id}",
        )
    conn = annotation_shard(annotation_id)
    try:
        try:
            page = list_namespace_accessions(
                conn,
                annotation_id,
                stored_ns,
                q=q,
                limit=limit,
                offset=offset,
            )
        except QueryError as exc:
            raise http_for_query_error(exc) from exc
        except ValueError as exc:
            raise HTTPException(
                status_code=422,
                detail=error_detail("invalid_offset", str(exc)),
            ) from exc
    finally:
        conn.close()
    if page is None:
        raise http_error(
            404,
            "namespace_not_found",
            f"Namespace {stored_ns!r} not found in {annotation_id}",
        )
    return page


@router.get("/annotations/{annotation_id}/genes")
def get_annotation_genes(
    annotation_id: str,
    names: str | None = Query(
        None,
        description=(
            "Comma-separated exact gene names (symbol∪alias). "
            "Wildcards are not supported. Mutually exclusive with curies."
        ),
    ),
    curies: str | None = Query(
        None,
        description=(
            "Comma-separated CURIEs (ns:accession) for exact xref resolve. "
            "Mutually exclusive with names."
        ),
    ),
    match: str = Query(MATCH_ANY, description=_MATCH_DESC),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
    cursor: str | None = Query(
        None, description="`next_cursor` of the previous page; use it instead of `offset`"),
):
    """Paginated gene cards for one annotation (optional names/curies filter)."""
    try:
        offset = resolve_offset(cursor, offset)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
    conn = annotation_shard(annotation_id)
    try:
        try:
            return list_annotation_genes(
                conn,
                annotation_id,
                names=names,
                curies=curies,
                match=match,
                limit=limit,
                offset=offset,
            )
        except QueryError as exc:
            raise http_for_query_error(exc) from exc
        except ValueError as exc:
            raise HTTPException(
                status_code=422,
                detail=error_detail("invalid_offset", str(exc)),
            ) from exc
    finally:
        conn.close()


@router.get("/annotations/{annotation_id}/genes/{local_id}/xrefs")
def get_annotation_gene_xrefs(
    annotation_id: str,
    local_id: int = Path(..., ge=0, le=2**63 - 1),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    offset: int = Query(0, ge=0, le=MAX_OFFSET),
    cursor: str | None = Query(
        None, description="`next_cursor` of the previous page; use it instead of `offset`"),
):
    """Paginated xref chips for one gene in this annotation's shard."""
    try:
        offset = resolve_offset(cursor, offset)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
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
                status_code=422,
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
    local_id: int = Path(..., ge=0, le=2**63 - 1),
):
    """One gene locus card from this annotation's gene corpus shard."""
    conn = annotation_shard(annotation_id)
    try:
        row = get_gene(conn, local_id, annotation_id)
    finally:
        conn.close()
    if row is None:
        raise http_error(
            404,
            "gene_not_found",
            f"Gene local_id {local_id} not found in {annotation_id}",
        )
    return row
