"""HTTP routes for annotrieve-xrefs (nginx strips /annotrieve-xrefs/api/v1)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Path, Query

from helpers.cursor import QueryError
from helpers.paths import MAX_ANNOTATION_ID_CHARS
from helpers.query import (
    ACCESSION_SORTS,
    DEFAULT_LIMIT,
    MAX_LIMIT,
    get_namespace,
    list_accessions,
    list_genes,
    list_namespaces,
)
from helpers.session import GeneCorpusUnavailable, open_reader
from helpers.shard_query import ShardUnavailable, get_gene_with_xrefs, open_shard

router = APIRouter()


def _conn():
    try:
        return open_reader()
    except GeneCorpusUnavailable as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc


@router.get("/health")
def health():
    return {"status": "ok"}


@router.get("/namespaces")
def get_namespaces():
    """Namespaces present in the published gene corpus."""
    return list_namespaces(_conn())


@router.get("/namespaces/{namespace:path}/accessions")
def get_namespace_accessions(
    namespace: str,
    sort: str = Query("n_loci"),
    sort_order: str = Query("desc"),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    cursor: str | None = Query(
        None,
        description="Opaque pagination cursor from a previous next/previous field",
    ),
):
    """Paginated accessions in one namespace."""
    if sort not in ACCESSION_SORTS:
        raise HTTPException(
            status_code=400,
            detail=f"sort must be one of: {', '.join(sorted(ACCESSION_SORTS))}",
        )
    if sort_order not in {"asc", "desc"}:
        raise HTTPException(status_code=400, detail="sort_order must be asc or desc")
    try:
        page = list_accessions(
            _conn(),
            namespace,
            sort=sort,
            sort_order=sort_order,
            cursor=cursor,
            limit=limit,
        )
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    if page is None:
        raise HTTPException(
            status_code=404, detail=f"Unknown namespace: {namespace}"
        )
    return page


@router.get("/namespaces/{namespace:path}")
def get_namespace_item(namespace: str):
    """One namespace summary."""
    row = get_namespace(_conn(), namespace)
    if row is None:
        raise HTTPException(
            status_code=404, detail=f"Unknown namespace: {namespace}"
        )
    return row


@router.get("/genes")
def get_genes(
    id: list[str] = Query(
        ...,
        description="CURIE(s) such as symbol:tp53 or GO:0008150 (repeatable)",
    ),
    taxid: int | None = Query(
        None,
        ge=1,
        description="Optional NCBI taxid (species or ancestor via annotation_lineage)",
    ),
    limit: int = Query(DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
    cursor: str | None = Query(
        None,
        description="Opaque pagination cursor from a previous next/previous field",
    ),
):
    """Gene rows matching any of the given CURIEs (keyset pagination)."""
    try:
        return list_genes(
            _conn(),
            id,
            taxid=taxid,
            cursor=cursor,
            limit=limit,
        )
    except QueryError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get("/annotations/{annotation_id}/genes/{local_id}")
def get_annotation_gene(
    annotation_id: str,
    local_id: int = Path(..., ge=0),
):
    """One gene locus card from this annotation's gene corpus shard."""
    if (
        not annotation_id
        or len(annotation_id) > MAX_ANNOTATION_ID_CHARS
        or "\0" in annotation_id
    ):
        raise HTTPException(status_code=404, detail="Gene shard not found")
    try:
        conn = open_shard(annotation_id)
    except GeneCorpusUnavailable as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except ShardUnavailable as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    try:
        row = get_gene_with_xrefs(conn, local_id)
    finally:
        conn.close()
    if row is None:
        raise HTTPException(
            status_code=404,
            detail=f"Gene local_id {local_id} not found in {annotation_id}",
        )
    return {"annotation_id": annotation_id, **row}
