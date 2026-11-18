"""Gene-name presence search (``/hits``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel, ConfigDict, Field

from helpers.concurrency import BusyError, hits_slot
from helpers.cursor import QueryError
from helpers.deps import corpus_conn
from helpers.errors import error_detail
from helpers.presence import presence_enabled, presence_reader
from helpers.pagination import MAX_NAMES_GET, MAX_NAMES_POST
from helpers.strict import reject_unknown_params
from services.hits import list_name_hits
from services.hits_formats import stream_hits_tsv
from services.presence_hits import list_name_hits_presence

router = APIRouter()

# Validation-style QueryError codes → HTTP 422
_VALIDATION_CODES = frozenset(
    {
        "names_required",
        "too_many_names",
        "too_many_taxids",
        "wildcard_not_supported",
        "invalid_taxid",
        "invalid_alias_mode",
        "invalid_source",
        "invalid_min_matched",
        "invalid_shape",
        "invalid_format",
        "presence_required",
        "invalid_terms",
        "unknown_term_field",
        "invalid_category",
        "option_needs_prefix",
        "invalid_prefix",
        "prefix_too_short",
        "prefix_too_broad",
    }
)
_FORMATS = ("json", "tsv")
_GET_PARAMS = ("names", "include_taxids", "exclude_taxids", "alias", "sources", "min_matched",
               "shape", "format")


class HitsSearchBody(BaseModel):
    """JSON body for ``POST /hits``.

    Cap is enforced in ``parse_name_list`` (casefold-deduped) so GET and POST
    share ``too_many_names`` — do not use list ``max_length`` here.
    """

    model_config = ConfigDict(extra="forbid")

    names: list[str] | None = Field(None, description="exact names or prefixes (seleno*)")
    terms: list[str | dict] | None = Field(
        None, description="terms as strings or objects {q, exclude, categories, split} "
                          "(options apply to prefix terms)")
    include_taxids: list[int] | None = None
    exclude_taxids: list[int] | None = None
    alias: str | None = None
    sources: list[str] | str | None = Field(
        None, description="Ensembl, GenBank and/or RefSeq; omit for all")
    min_matched: int | None = Field(
        None, description="keep annotations that match at least this many terms (symbol hits)")
    shape: str | None = Field(None, description="compact (default) | records")
    format: str | None = Field(None, description="json (default) | tsv (implies records)")


def _busy_response() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={
            "detail": error_detail(
                "hits_busy",
                "too many concurrent search requests; retry shortly",
            )
        },
        headers={"Retry-After": "1"},
    )


def _http_for_query_error(exc: QueryError) -> HTTPException:
    status = 422 if exc.code in _VALIDATION_CODES else 400
    return HTTPException(status_code=status, detail=exc.as_detail())


def _run(*, names, terms=None, include_taxids, exclude_taxids, alias, sources, min_matched, shape, fmt, max_names):
    """One search for GET and POST; returns a JSON-able dict or a streaming TSV response."""
    fmt = (fmt or "json").strip().lower()
    if fmt not in _FORMATS:
        raise QueryError("format must be json or tsv", code="invalid_format", value=fmt)
    if not presence_enabled():
        if sources or min_matched or (shape or "compact") != "compact" or fmt != "json":
            raise QueryError("sources, min_matched, shape and format need the presence index",
                             code="presence_required")
        if terms:
            raise QueryError("terms need the presence index", code="presence_required")
        return list_name_hits(corpus_conn(), names, include_taxids=include_taxids,
                              exclude_taxids=exclude_taxids, max_names=max_names)
    reader = presence_reader()
    result = list_name_hits_presence(
        reader, names, terms=terms, include_taxids=include_taxids, exclude_taxids=exclude_taxids, alias=alias,
        max_names=max_names, sources=sources, min_matched=min_matched,
        shape="records" if fmt == "tsv" else shape)
    if fmt == "tsv":
        return StreamingResponse(
            stream_hits_tsv(result), media_type="text/tab-separated-values; charset=utf-8",
            headers={"X-Data-Version": str(result["data_version"])})
    return result


@router.get("/hits")
def get_hits(
    request: Request,
    names: str = Query(
        ...,
        description=(
            "Comma-separated gene names (symbol∪alias); a trailing * makes a prefix term "
            "(seleno*, at least 3 characters, expanded server-side and echoed in `terms`)."
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
            "off | suggest | include. suggest (default) counts symbol hits as matches and "
            "lists alias-only annotations as `alias_only`; include restores the union. "
            "Only applies when the presence index is enabled (PRESENCE_DIR)."
        ),
    ),
    sources: str | None = Query(
        None, description="Comma-separated Ensembl, GenBank, RefSeq; omit for all"),
    min_matched: int | None = Query(
        None, ge=1, description="Keep annotations that match at least this many terms"),
    shape: str = Query(
        "compact", description="compact (positional fields/rows) | records (one object per annotation)"),
    format: str = Query(
        "json", description="json | tsv (a table with # provenance lines; implies records)"),
):
    """
    Presence search by exact gene symbol/synonym names.

    ``shape=compact``: ``terms``, ``fields``/``rows`` with ``matched`` (term indices).
    ``shape=records``: one object per annotation with ``n_matched``, ``n_terms``, ``missing`` and
    ``download_url``. Both carry ``facets`` (annotation counts per source, before the ``sources``
    filter is applied).
    """
    reject_unknown_params(request, _GET_PARAMS)
    try:
        with hits_slot():
            return _run(names=names, include_taxids=include_taxids, exclude_taxids=exclude_taxids,
                        alias=alias, sources=sources, min_matched=min_matched, shape=shape,
                        fmt=format, max_names=MAX_NAMES_GET)
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc


@router.post("/hits")
def post_hits(body: HitsSearchBody):
    """Presence search for a larger name list (JSON body). Prefer GET when ≤200."""
    try:
        with hits_slot():
            return _run(names=body.names, terms=body.terms, include_taxids=body.include_taxids,
                        exclude_taxids=body.exclude_taxids, alias=body.alias, sources=body.sources,
                        min_matched=body.min_matched, shape=body.shape, fmt=body.format,
                        max_names=MAX_NAMES_POST)
    except BusyError:
        return _busy_response()
    except QueryError as exc:
        raise _http_for_query_error(exc) from exc
