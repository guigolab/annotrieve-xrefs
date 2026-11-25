"""Gene length per annotation, by clade (``/distributions``)."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Path, Query, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict

from helpers.concurrency import BusyError, acquire_loci_slot, release_loci_slot
from helpers.cursor import QueryError
from helpers.errors import error_detail, http_error, http_for_query_error
from helpers.pagination import MAX_NAMES_POST
from helpers.presence import (
    annotation_info,
    annotation_paths,
    open_request_connections,
    presence_enabled,
    presence_reader,
)
from helpers.strict import reject_unknown_params
from services.distributions import (
    BINS,
    POINT_FIELDS,
    SUMMARY_FIELDS,
    column_masks,
    names_of,
    parse_columns,
    points,
    summarise,
)
from services.presence_hits import SOURCE_LABELS
from services.presence_loci import build_query, estimate
from services.terms import term_echo
from settings import settings

router = APIRouter()

_POINTS_PARAMS = ("name", "alias", "sources", "annotation_ids")


class DistributionColumn(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str | None = None
    taxids: list[int] | None = None


class DistributionsBody(BaseModel):
    """Names (or prefix terms), the clades to group by, and the usual filters."""

    model_config = ConfigDict(extra="forbid")

    names: list[str] | None = None
    terms: list[str | dict] | None = None
    columns: list[DistributionColumn] | None = None
    include_taxids: list[int] | None = None
    exclude_taxids: list[int] | None = None
    alias: str | None = "suggest"
    annotation_ids: list[str] | None = None
    sources: list[str] | str | None = None
    dry_run: bool = False


def _busy() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={"detail": error_detail("loci_busy", "too many concurrent requests; retry shortly")},
        headers={"Retry-After": "1"},
    )


def _head(q, data_version: str) -> dict:
    return {"data_version": data_version, "terms": [term_echo(t) for t in q.terms],
            "alias": q.mode, "include_taxids": q.include, "exclude_taxids": q.exclude,
            "source_filter": [SOURCE_LABELS[c] for c in q.sources]}


def _run_summary(body: DistributionsBody):
    if not presence_enabled():
        raise http_error(503, "distributions_unavailable", "The presence/loci index is not enabled")
    columns = parse_columns(body.columns)
    reader = presence_reader()
    try:
        acquire_loci_slot()
    except BusyError:
        return _busy()
    pres = loci = None
    try:
        pres, loci, data_version = open_request_connections()
        q = build_query(pres, reader.all_mask, data_version, body.names, terms=body.terms,
                        include_taxids=body.include_taxids, exclude_taxids=body.exclude_taxids,
                        alias=body.alias, annotation_ids=body.annotation_ids, sources=body.sources,
                        max_names=MAX_NAMES_POST)
        names = names_of(q)
        est = estimate(pres, loci, q)
        head = {**_head(q, data_version), "bins": BINS,
                "names": [{"name": n["name"], "term": n["term"]} for n in names], **est}
        col_sets = column_masks(pres, reader.all_mask, columns, q.scope)
        head["columns"] = [{"id": c["id"], "taxids": c["taxids"], "n_annotations": len(s)}
                           for c, s in zip(columns, col_sets)]
        if body.dry_run:
            return head
        cap = int(settings.DISTRIBUTIONS_MAX_LOCI)
        if cap > 0 and est["n_loci_upper_bound"] > cap:
            raise http_error(413, "distributions_too_large",
                             "this request would read more genes than allowed; ask for fewer names "
                             "or narrow it with include_taxids / sources",
                             n_loci_upper_bound=est["n_loci_upper_bound"], max=cap)
        rows = summarise(loci, q, names, col_sets)
        return {**head, "fields": list(SUMMARY_FIELDS), "total": len(rows), "rows": rows}
    finally:
        for conn in (pres, loci):
            if conn is not None:
                conn.close()
        release_loci_slot()


@router.post("/distributions", summary="Gene length per annotation, summarised by clade")
def post_distributions(body: DistributionsBody):
    """Length of each gene (end - start + 1) across annotations, grouped into the clades you give.

    ``columns`` is a list of ``{id, taxids}``; each column holds the annotations under any of its
    taxids (no columns: one column with everything). The answer is a compact table like `/hits`:
    ``fields`` names the positions of each row, and a row is one (name, column) cell with the
    number of annotations that have the gene, how many of them have several copies, the exact
    mean / min / quartiles / max of the per-annotation length, and a 40-bucket histogram over the
    fixed log scale in ``bins``. Cells without the gene are left out; ``columns[].n_annotations``
    is the denominator.

    When an annotation has several matching genes, its value is their mean length and it is
    counted in ``n_multi_copy``. Genes that only carry the name as an alias are used only when
    the annotation has no symbol gene, and never with ``alias=off``.
    """
    try:
        return _run_summary(body)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _run_points(taxid: int, name: str, alias, sources, annotation_ids):
    if not presence_enabled():
        raise http_error(503, "distributions_unavailable", "The presence/loci index is not enabled")
    reader = presence_reader()
    try:
        acquire_loci_slot()
    except BusyError:
        return _busy()
    pres = loci = None
    try:
        pres, loci, data_version = open_request_connections()
        q = build_query(pres, reader.all_mask, data_version, name, include_taxids=[taxid],
                        alias=alias, annotation_ids=annotation_ids, sources=sources, max_names=1)
        term = q.terms[0]
        if term["kind"] != "exact":
            raise QueryError("name must be one exact gene name, not a prefix",
                             code="exact_name_required", value=name)
        rows = []
        if term["name_ids"]:
            rows = points(loci, q, term["name_ids"][0], annotation_info(reader),
                          annotation_paths(reader))
        return {**_head(q, data_version), "taxid": taxid, "name": term["name"],
                "n_annotations": len(q.scope), "fields": list(POINT_FIELDS),
                "total": len(rows), "rows": rows}
    finally:
        for conn in (pres, loci):
            if conn is not None:
                conn.close()
        release_loci_slot()


@router.get("/distributions/{taxid}/points", summary="Gene length of every annotation of one clade")
def get_distribution_points(
    request: Request,
    taxid: int = Path(..., ge=1, description="NCBI taxid of the clade"),
    name: str = Query(..., description="One exact gene name"),
    alias: str = Query("suggest", description="off | suggest | include; off drops alias-only genes"),
    sources: str | None = Query(None, description="Comma-separated Ensembl, GenBank, RefSeq"),
    annotation_ids: str | None = Query(None, description="Restrict to these annotation ids"),
):
    """One row per annotation under ``taxid`` that has the gene: ``[annotation_id, taxid, length,
    copies, via]``. ``length`` is the mean when there are several copies (``copies`` > 1); ``via`` is
    0 for symbol genes and 1 when only an alias matched. Same envelope as ``POST /distributions``,
    so a client can draw a clade's cell and its dots from the same fields. Organism, assembly and
    source are not repeated: join them on ``annotation_id`` from ``/hits``."""
    reject_unknown_params(request, _POINTS_PARAMS)
    try:
        return _run_points(taxid, name, alias, sources, annotation_ids)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
