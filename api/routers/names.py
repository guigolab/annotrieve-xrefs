"""Gene-name prefix discovery (``/names``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.presence import presence_enabled, presence_reader
from services.names import list_names_by_prefix
from services.presence_hits import list_names_by_prefix_presence

router = APIRouter()

_VALIDATION_CODES = frozenset(
    {
        "prefix_too_short",
        "prefix_too_broad",
        "invalid_prefix",
        "too_many_taxids",
        "invalid_taxid",
        "invalid_alias_mode",
    }
)


def _http_for_query_error(exc: QueryError) -> HTTPException:
    status = 422 if exc.code in _VALIDATION_CODES else 400
    return HTTPException(status_code=status, detail=exc.as_detail())


@router.get("/names")
def get_names(
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
):
    """
    Prefix discovery for building an exact gene-name set.

    Returns categorized accessions with taxon-scoped annotation counts.
    Use ``GET|POST /hits`` with the selected exact names to find annotations.
    """
    try:
        if presence_enabled():
            return list_names_by_prefix_presence(
                presence_reader(), prefix, include_taxids=include_taxids,
                exclude_taxids=exclude_taxids, alias=alias)
        return list_names_by_prefix(
            corpus_conn(),
            prefix,
            include_taxids=include_taxids,
            exclude_taxids=exclude_taxids,
        )
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc
