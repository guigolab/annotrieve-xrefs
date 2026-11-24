"""Gene-name prefix discovery (``/names``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Request

from helpers.cursor import QueryError
from helpers.errors import http_for_query_error
from helpers.strict import reject_unknown_params
from helpers.presence import open_request_connections, presence_reader
from services.presence_hits import list_names_by_prefix_presence

router = APIRouter()



def _open_loci():
    """loci.sqlite connection for ``n_loci``, or None when that index is not published."""
    try:
        pres, loci, _v = open_request_connections()
    except HTTPException:
        return None
    pres.close()
    return loci


_GET_PARAMS = ("prefix", "include_taxids", "exclude_taxids", "alias", "sources")


@router.get("/names")
def get_names(
    request: Request,
    prefix: str = Query(
        ...,
        description=(
            "Gene-name prefix stem (symbol∪alias); no wildcard character; "
            "at least 3 characters"
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
    alias: str = Query(
        "suggest",
        description=(
            "off | suggest | include. n_annotations counts symbol annotations; "
            "n_alias_extra counts annotations reachable only through an alias. "
            "Only applies when the presence index is enabled (PRESENCE_DIR)."
        ),
    ),
    sources: str | None = Query(
        None, description="Comma-separated Ensembl, GenBank, RefSeq; counts only those annotations"),
):
    """
    Prefix discovery for building an exact gene-name set.

    Returns categorized accessions with taxon-scoped annotation counts.
    Use ``GET|POST /hits`` with the selected exact names to find annotations.
    """
    reject_unknown_params(request, _GET_PARAMS)
    try:
        loci = _open_loci()
        try:
            return list_names_by_prefix_presence(
                presence_reader(), prefix, include_taxids=include_taxids,
                exclude_taxids=exclude_taxids, alias=alias, sources=sources, loci=loci)
        finally:
            if loci is not None:
                loci.close()
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
