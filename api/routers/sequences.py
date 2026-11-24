"""FASTA export (``/sequences``), sequence estimates and plans, jobs and cached results.

See API-DESIGN.md section 4.4 for the contract. ``result_id`` identifies a request (normalised
body + data version): equal requests reuse one result on disk.
"""
from __future__ import annotations

import logging
import shutil
import time
from typing import Any

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, ConfigDict, Field

from helpers.client_ip import client_ip, owner_key
from helpers.concurrency import BusyError, sequences_slot
from helpers.cursor import QueryError
from helpers.errors import error_detail, http_for_query_error
from services.seq_export import count_sequences, run_export
from services.seq_request import SeqRequest, parse_request
from services.seq_resolve import GeneResolver, SeqError
from services.seq_results import (
    Refused,
    admit,
    purge_expired,
    read_status,
    result_dir,
    valid_id,
    write_status,
)
from services.seq_source import Source
from settings import settings

router = APIRouter()

_ERR = {"description": "Error: `detail` is `{code, message, ...}`."}
_FASTA_HEADERS = {
    "X-Sequence-Count": {"description": "Sequences in the FASTA", "schema": {"type": "integer"}},
    "X-Skipped-Count": {"description": "Genes skipped (see the manifest)", "schema": {"type": "integer"}},
    "X-Gene-Count": {"description": "Genes processed", "schema": {"type": "integer"}},
    "X-Result-Id": {"description": "Identifier of this result", "schema": {"type": "string"}},
    "X-Data-Version": {"description": "Index version used", "schema": {"type": "string"}},
    "X-Manifest-Url": {"description": "Per-gene outcome table (tsv)", "schema": {"type": "string"}},
}
SEQUENCES_RESPONSES = {
    200: {"description": "FASTA (`text/x-fasta`). Headers: result id, counts, manifest link.",
          "headers": _FASTA_HEADERS,
          "content": {"text/x-fasta": {"schema": {"type": "string"}}}},
    404: {**_ERR, "description": "`no_genes_matched`, or a pinned/unknown locus (loci mode, on_error=fail)"},
    409: {**_ERR, "description": "`export_running`: an identical request is being produced; retry shortly"},
    413: {**_ERR, "description": "`too_many_genes` or `too_many_sequences` (over the limit of one request)"},
    422: {**_ERR, "description": "Invalid request (`validation_error`, `exactly_one_source`, "
                                 "`unknown_rank_key`, `invalid_locus_ref`, ...)"},
    429: {**_ERR, "description": "`too_many_active_exports` or `export_rate_limited` for this client "
                                 "IP (see `Retry-After`)"},
    503: {**_ERR, "description": "`sequences_busy`, `sequences_capacity`, `annotations_unavailable`, "
                                 "`locus_ref_index_unavailable`, `gffread_unavailable`"},
    504: {**_ERR, "description": "`sequences_timeout`"},
}


logger = logging.getLogger(__name__)


class SequenceEstimate(BaseModel):
    """What a request would cost, from the index alone (nothing is extracted)."""

    data_version: str
    result_id: str = Field(description="Id the finished result would have")
    n_genes: int | None = Field(None, description="null when over the limit (see n_genes_at_least)")
    n_genes_at_least: int
    complete: bool = Field(description="false when counting stopped at max_genes")
    est_bp_upper_bound: int = Field(description="Sum of gene lengths: bound of the bases to fetch")
    max_genes: int
    over_limit: bool
    max_sequences: int = Field(description="most sequences one request may return")
    cached: bool = Field(description="a finished result for this request exists")
    terms: list[dict[str, Any]] = Field(description="query mode: each name and whether it was found")


class PlanTranscript(BaseModel):
    model_config = ConfigDict(extra="allow")

    transcript_id: str
    selected: bool
    coding: bool
    cds_len: int
    spliced_len: int
    span: int
    n_exons: int
    utr_len: int
    biotype: str | None = None
    tags: list[str] = []


class PlanGene(BaseModel):
    model_config = ConfigDict(extra="allow")

    locus_ref: str
    status: str = Field(description="ok | skipped")
    reason: str = Field("", description="skip reason code when status is skipped")
    message: str = ""
    n_candidates: int = 0
    transcripts: list[PlanTranscript] = []


class PlanLocus(BaseModel):
    locus_ref: str
    transcript_id: str


class SequencePlan(BaseModel):
    """The transcripts a request would extract, with the metrics behind the choice."""

    data_version: str
    action: str
    n_genes: int
    n_skipped: int
    genes: list[PlanGene]
    loci: list[PlanLocus] = Field(description="the selection, ready to post to /sequences as `loci`")


class SequencesBody(BaseModel):
    """Exactly one of ``query``, ``loci`` or ``items`` selects the genes."""

    model_config = ConfigDict(extra="forbid")     # a typo'd field is an error, not silently ignored

    action: str | None = Field(None, description="cds | transcript | protein")
    query: dict[str, Any] | None = Field(
        None, description="names, include_taxids, exclude_taxids, alias, annotation_ids, "
                          "case_sensitive (same meaning as /loci)")
    loci: list[Any] | None = Field(
        None, description="'annotation_id:local_id' strings or {locus_ref, transcript_id} objects")
    items: list[Any] | None = Field(
        None, description="legacy: {annotation_id, local_id, transcript_id} objects")
    transcript: str | dict[str, Any] | None = Field(
        None, description="preset (longest, longest_spliced, longest_cds, shortest_cds, "
                          "most_exons, fewest_exons, all) or {rank_by, require, keep}")
    buffer_bp: int = 0
    include_exon_coords: bool = False
    on_error: str | None = Field(None, description="skip (default for query) | fail (default for loci)")
    max_genes: int | None = None


def _http(exc: SeqError) -> HTTPException:
    headers = {"Retry-After": str(exc.retry_after)} if exc.retry_after is not None else None
    return HTTPException(status_code=exc.status,
                         detail=error_detail(exc.code, str(exc), **exc.extra), headers=headers)


def _busy() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={"detail": error_detail("sequences_busy",
                                        "too many concurrent sequences requests; retry shortly")},
        headers={"Retry-After": "1"})


def _parse(body: SequencesBody, *, default_action: str | None = None) -> SeqRequest:
    data = body.model_dump(exclude_none=True)
    if default_action and "action" not in data:
        data["action"] = default_action
    try:
        return parse_request(data)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _urls(request: Request, rid: str) -> dict[str, str]:
    base = (settings.PUBLIC_BASE_URL or str(request.base_url)).rstrip("/")
    return {"result_url": f"{base}/sequences/results/{rid}",
            "manifest_url": f"{base}/sequences/results/{rid}/manifest.tsv"}


def _open_source(req: SeqRequest) -> Source:
    try:
        return Source(req)
    except SeqError as exc:
        raise _http(exc) from exc
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _count(src: Source, req: SeqRequest):
    try:
        return src.count(req.max_genes)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _fasta_response(request: Request, rid: str, status: dict) -> FileResponse:
    d = result_dir(rid)
    return FileResponse(
        d / "result.fa", media_type="text/x-fasta", filename=f"annotrieve-{rid[:8]}.fa",
        headers={"X-Sequence-Count": str(status.get("n_sequences", 0)),
                 "X-Skipped-Count": str(status.get("n_skipped", 0)),
                 "X-Gene-Count": str(status.get("n_genes", 0)),
                 "X-Result-Id": rid, "X-Data-Version": str(status.get("data_version", "")),
                 "X-Manifest-Url": _urls(request, rid)["manifest_url"],
                 "Cache-Control": "no-store"})


def _summary(result, data_version: str) -> dict:
    return {"n_sequences": result.n_sequences, "n_genes": result.n_genes,
            "n_skipped": result.n_skipped, "reasons": result.reasons,
            "data_version": data_version, "error": None}


@router.post("/sequences", responses=SEQUENCES_RESPONSES, summary="FASTA for a query or a list of loci")
def post_sequences(body: SequencesBody, request: Request):
    """Extract CDS, transcript or protein FASTA for the genes of a query or an explicit list."""
    req = _parse(body)
    src = _open_source(req)
    try:
        rid = req.result_id(src.data_version)
        purge_expired()
        status = read_status(rid)
        if status and status["state"] == "done" and (result_dir(rid) / "result.fa").is_file():
            return _fasta_response(request, rid, status)
        if status and status["state"] in ("queued", "running"):
            raise HTTPException(
                status_code=409, headers={"Retry-After": "5"},
                detail=error_detail("export_running", "an identical request is being produced; "
                                    "retry in a few seconds (it will then be served from the store)"))

        counted = _count(src, req)
        if not counted.complete:
            raise HTTPException(status_code=413, detail=error_detail(
                "too_many_genes", f"the request matches more than {req.max_genes} genes; narrow "
                "it or lower max_genes", n_genes_at_least=counted.n_genes, max_genes=req.max_genes))
        if counted.n_genes == 0:
            raise HTTPException(status_code=404, detail=error_detail(
                "no_genes_matched", "no gene matches the request", terms=src.terms()))
        if req.can_give_several_per_gene:
            # One gene can give several sequences: count them (local GFF reads, no network) and
            # refuse before the first NCBI call.
            limit = int(settings.SEQUENCES_MAX_SEQUENCES)
            n_seq = count_sequences(req, src.items(), limit)
            if n_seq > limit:
                raise HTTPException(status_code=413, detail=error_detail(
                    "too_many_sequences", f"the request gives more than {limit} sequences; narrow it "
                    "or choose one transcript per gene", n_sequences_at_least=n_seq,
                    max_sequences=limit))
        data_version = src.data_version

        owner = owner_key(client_ip(request))
        try:
            admit(owner)
        except Refused as exc:
            raise HTTPException(
                status_code=exc.status, detail=error_detail(exc.code, str(exc), **exc.extra),
                headers={"Retry-After": str(exc.retry_after)}) from exc

        try:
            with sequences_slot():
                write_status(rid, state="running", done=0, total=counted.n_genes,
                             created=time.time(), owner=owner)
                t0 = time.monotonic()
                try:
                    result = run_export(req, src.items(), out_dir=result_dir(rid),
                                        total=counted.n_genes)
                except SeqError as exc:
                    shutil.rmtree(result_dir(rid), ignore_errors=True)
                    raise _http(exc) from exc
                except QueryError as exc:
                    shutil.rmtree(result_dir(rid), ignore_errors=True)
                    raise http_for_query_error(exc) from exc
                except BaseException:
                    # Never leave a "running" status behind: equal requests would be told
                    # to wait for an export that does not exist.
                    shutil.rmtree(result_dir(rid), ignore_errors=True)
                    raise
                status = write_status(rid, state="done", **_summary(result, data_version))
                logger.info("POST /sequences genes=%d seqs=%d skipped=%d wall_s=%.2f",
                            result.n_genes, result.n_sequences, result.n_skipped,
                            time.monotonic() - t0)
        except BusyError:
            return _busy()
        return _fasta_response(request, rid, status)
    finally:
        src.close()


@router.post("/sequence-estimates", response_model=SequenceEstimate,
             summary="Estimate a sequence request", responses={413: _ERR, 422: _ERR, 503: _ERR})
def create_sequence_estimate(body: SequencesBody, request: Request):
    """Cheap preview: gene count and an upper bound of the bases to fetch. Reads only the index."""
    req = _parse(body, default_action="cds")
    src = _open_source(req)
    try:
        counted = _count(src, req)
        rid = req.result_id(src.data_version)
        status = read_status(rid)
        return {
            "data_version": src.data_version, "result_id": rid,
            "n_genes": counted.n_genes if counted.complete else None,
            "n_genes_at_least": counted.n_genes, "complete": counted.complete,
            "est_bp_upper_bound": counted.bp_upper_bound, "max_genes": req.max_genes,
            "over_limit": not counted.complete,
            "max_sequences": int(settings.SEQUENCES_MAX_SEQUENCES),
            "cached": bool(status and status["state"] == "done"),
            "terms": src.terms(),
        }
    finally:
        src.close()


@router.post("/sequence-plans", response_model=SequencePlan,
             summary="Plan the transcripts of a sequence request",
             responses={413: _ERR, 422: _ERR, 503: _ERR})
def create_sequence_plan(body: SequencesBody):
    """Which transcripts the request would extract, with the metrics behind the choice.

    Reads the GFF slices (no network). A client can edit the result and post the selected
    ``loci`` back to ``/sequences``.
    """
    req = _parse(body)
    src = _open_source(req)
    try:
        cap = min(req.max_genes, int(settings.SEQUENCES_PLAN_MAX_GENES))
        counted = _count_cap(src, cap)
        if not counted.complete:
            raise HTTPException(status_code=413, detail=error_detail(
                "plan_too_large", f"a plan covers at most {cap} genes; narrow the request",
                n_genes_at_least=counted.n_genes, max_genes=cap))
        resolver = GeneResolver(req.action.value, req.rule, req.buffer_bp)
        genes, loci_out = [], []
        for ref, gene, reason in _plan_items(src, resolver):
            if gene is None:
                genes.append({"locus_ref": ref, "status": "skipped", "reason": reason,
                              "message": reason.replace("_", " "), "transcripts": []})
                continue
            try:
                res = resolver.resolve(gene)
            except SeqError as exc:
                raise _http(exc) from exc
            chosen = {c.transcript_id for c in res.chosen}
            genes.append({
                "locus_ref": ref, "annotation_id": gene.annotation_id, "symbol": gene.symbol,
                "organism_name": gene.organism, "assembly_accession": gene.accession,
                "seqid": gene.seqid, "start": gene.start, "end": gene.end, "strand": gene.strand,
                "status": "ok" if res.ok else "skipped", "reason": res.reason,
                "message": res.message, "rule": res.rule, "n_candidates": res.n_candidates,
                "transcripts": [{"transcript_id": c.transcript_id,
                                 "selected": c.transcript_id in chosen, **c.metrics()}
                                for c in res.candidates]})
            loci_out += [{"locus_ref": ref, "transcript_id": t} for t in sorted(chosen)]
        return {"data_version": src.data_version, "action": req.action.value,
                "n_genes": len(genes), "n_skipped": sum(g["status"] == "skipped" for g in genes),
                "genes": genes, "loci": loci_out}
    finally:
        src.close()


def _plan_items(src: Source, resolver: GeneResolver):
    """The source's genes; releases the resolver's open GFF files when the loop ends."""
    try:
        yield from src.items()
    finally:
        resolver.close()


def _count_cap(src: Source, cap: int):
    try:
        return src.count(cap)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _ready(result_id: str) -> tuple[dict, Any]:
    status = read_status(result_id) if valid_id(result_id) else None
    if status is None:
        raise HTTPException(status_code=404, detail=error_detail(
            "job_not_found", "no such result (results expire after a while)"))
    if status["state"] == "failed":
        err = status.get("error") or {}
        raise HTTPException(status_code=int(err.get("status", 500)), detail=error_detail(
            err.get("code", "job_failed"), err.get("message", "the job failed")))
    if status["state"] != "done":
        raise HTTPException(status_code=409, detail=error_detail(
            "result_not_ready", "the job has not finished", state=status["state"],
            done=status.get("done"), total=status.get("total")))
    d = result_dir(result_id)
    if not (d / "result.fa").is_file():
        raise HTTPException(status_code=404, detail=error_detail(
            "job_not_found", "this result has expired; repeat the request"))
    return status, d


@router.get("/sequences/results/{result_id}", summary="Download a finished FASTA",
            responses={200: SEQUENCES_RESPONSES[200], 404: _ERR, 409: {**_ERR, "description": "`result_not_ready`"}})
def get_sequence_result(result_id: str, request: Request):
    status, _ = _ready(result_id)
    return _fasta_response(request, result_id, status)


@router.get("/sequences/results/{result_id}/manifest.tsv", summary="Per-gene outcome table",
            responses={200: {"description": "tsv", "content": {"text/tab-separated-values": {"schema": {"type": "string"}}}}, 404: _ERR, 409: _ERR})
def get_sequence_manifest(result_id: str):
    _status, d = _ready(result_id)
    return FileResponse(d / "manifest.tsv", media_type="text/tab-separated-values",
                        filename=f"annotrieve-{result_id[:8]}-manifest.tsv",
                        headers={"Cache-Control": "no-store"})
