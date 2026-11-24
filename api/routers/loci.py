"""Gene coordinates (``/loci``): where the matching genes are."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel, ConfigDict

from helpers.concurrency import BusyError, acquire_loci_slot, release_loci_slot
from helpers.cursor import QueryError
from helpers.errors import error_detail, http_error, http_for_query_error
from helpers.pagination import MAX_NAMES_GET, MAX_NAMES_POST, clamp_limit
from helpers.strict import reject_unknown_params
from helpers.presence import (
    annotation_info,
    open_request_connections,
    presence_enabled,
    presence_reader,
)
from services.loci_formats import stream_ndjson, stream_tsv
from services.presence_hits import SOURCE_LABELS
from services.presence_loci import (
    build_query,
    decode_cursor,
    encode_cursor,
    estimate,
    iter_loci,
)
from services.terms import term_echo
from settings import settings

router = APIRouter()

FORMATS = ("json", "tsv", "ndjson")
HARD_ROW_CEILING = 20_000_000
_GET_PARAMS = ("names", "include_taxids", "exclude_taxids", "alias", "annotation_ids",
               "sources", "case_sensitive", "flank_bp", "format", "dry_run", "limit", "cursor")
_MEDIA = {"tsv": "text/tab-separated-values", "ndjson": "application/x-ndjson"}


class LociBody(BaseModel):
    model_config = ConfigDict(extra="forbid")

    names: list[str] | None = None
    terms: list[str | dict] | None = None
    include_taxids: list[int] | None = None
    exclude_taxids: list[int] | None = None
    alias: str | None = "suggest"
    annotation_ids: list[str] | None = None
    sources: list[str] | str | None = None
    case_sensitive: bool = False
    flank_bp: int = 0
    format: str = "json"
    dry_run: bool = False
    limit: int | None = None
    cursor: str | None = None


def _busy() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={"detail": error_detail("loci_busy", "too many concurrent loci requests; retry shortly")},
        headers={"Retry-After": "1"},
    )


def _run(names, *, terms=None, include_taxids, exclude_taxids, alias, annotation_ids, sources,
         case_sensitive, flank_bp, fmt, dry_run, limit, cursor, max_names):
    if not presence_enabled():
        raise http_error(503, "loci_unavailable", "The presence/loci index is not enabled")
    fmt = (fmt or "json").lower()
    if fmt not in FORMATS:
        raise QueryError(f"format must be one of {', '.join(FORMATS)}", code="invalid_format",
                         value=fmt)
    reader = presence_reader()
    try:
        acquire_loci_slot()
    except BusyError:
        return _busy()
    pres = loci = None
    handed_off = False
    try:
        pres, loci, data_version = open_request_connections()
        q = build_query(pres, reader.all_mask, data_version, names, terms=terms,
                        include_taxids=include_taxids, exclude_taxids=exclude_taxids, alias=alias,
                        annotation_ids=annotation_ids, sources=sources, case_sensitive=case_sensitive,
                        flank_bp=flank_bp, max_names=max_names)
        est = estimate(pres, loci, q)
        info = annotation_info(reader)
        base = {"terms": [term_echo(t) for t in q.terms], "alias": q.mode, "include_taxids": q.include,
                "exclude_taxids": q.exclude,
                "source_filter": [SOURCE_LABELS[c] for c in q.sources], "data_version": data_version, **est}
        if dry_run:
            return base
        if fmt == "json":
            page = clamp_limit(limit if limit is not None else 100)
            after = decode_cursor(q, cursor)
            rows = []
            cur_ord, cur_count = (after if after else (None, 0))
            for row in iter_loci(pres, loci, q, after, info):
                if row["ord"] != cur_ord:
                    cur_ord, cur_count = row["ord"], 0
                if len(rows) >= page:
                    base["next_cursor"] = encode_cursor(q, cur_ord, cur_count)
                    break
                cur_count += 1
                rows.append({k: v for k, v in row.items() if k != "ord"})
            base.setdefault("next_cursor", None)
            return {**base, "rows": rows}

        cap = settings.LOCI_MAX_ROWS if limit is None else max(1, min(int(limit), HARD_ROW_CEILING))
        if limit is None and not q.scoped and est["n_loci_upper_bound"] > cap:
            raise http_error(413, "export_too_large",
                             "this query may return more rows than the export limit; narrow it "
                             "or raise limit=", n_loci_upper_bound=est["n_loci_upper_bound"],
                             max=cap)
        writer = {"tsv": stream_tsv, "ndjson": stream_ndjson}[fmt]
        rows_iter = iter_loci(pres, loci, q, None, info)

        def body():
            try:
                yield from writer(q, rows_iter, cap)
            finally:
                pres.close()
                loci.close()
                release_loci_slot()

        handed_off = True
        return StreamingResponse(body(), media_type=_MEDIA[fmt],
                                 headers={"X-Accel-Buffering": "no", "X-Data-Version": data_version})
    finally:
        if not handed_off:
            for conn in (pres, loci):
                if conn is not None:
                    conn.close()
            release_loci_slot()


def _wrap(**kwargs):
    try:
        return _run(**kwargs)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


@router.get("/loci")
def get_loci(
    request: Request,
    names: str = Query(..., description="Comma-separated gene names; name* is a prefix term (3+ characters, expanded server-side, echoed in `terms`)"),
    include_taxids: str | None = Query(None, description="Comma-separated NCBI taxids (lineage OR)"),
    exclude_taxids: str | None = Query(None, description="Comma-separated NCBI taxids to exclude"),
    alias: str = Query("suggest", description="off | suggest | include; off drops alias-only genes"),
    annotation_ids: str | None = Query(None, description="Restrict to these annotation ids"),
    sources: str | None = Query(None, description="Comma-separated Ensembl, GenBank, RefSeq; omit for all"),
    case_sensitive: bool = Query(False, description="Exact match on the gene's own spelling"),
    flank_bp: int = Query(0, description="Extra bases on each side of the gff_url slice"),
    format: str = Query("json", description="json | tsv | ndjson"),
    dry_run: bool = Query(False, description="Return counts only"),
    limit: int | None = Query(None, description="JSON: page size. Other formats: stop after this many rows (omit to use the export cap, which refuses oversized unscoped queries with 413)"),
    cursor: str | None = Query(None, description="JSON pagination cursor"),
):
    """Gene coordinates for the matching genes, one row per gene, deterministic order."""
    reject_unknown_params(request, _GET_PARAMS)
    return _wrap(names=names, include_taxids=include_taxids, exclude_taxids=exclude_taxids,
                 alias=alias, annotation_ids=annotation_ids, sources=sources,
                 case_sensitive=case_sensitive, flank_bp=flank_bp, fmt=format, dry_run=dry_run, limit=limit, cursor=cursor,
                 max_names=MAX_NAMES_GET)


@router.post("/loci")
def post_loci(body: LociBody):
    """Same as GET /loci for longer name lists."""
    return _wrap(names=body.names, terms=body.terms, include_taxids=body.include_taxids,
                 exclude_taxids=body.exclude_taxids, alias=body.alias,
                 annotation_ids=body.annotation_ids, sources=body.sources,
                 case_sensitive=body.case_sensitive,
                 flank_bp=body.flank_bp, fmt=body.format, dry_run=body.dry_run, limit=body.limit,
                 cursor=body.cursor, max_names=MAX_NAMES_POST)
