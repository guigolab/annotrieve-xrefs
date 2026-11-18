"""FASTA export (``/sequences``), its estimate/plan, jobs and cached results.

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
from pydantic import BaseModel, Field

from helpers.concurrency import BusyError, sequences_slot
from helpers.cursor import QueryError
from helpers.errors import error_detail
from services.seq_export import run_export
from services.seq_request import SeqRequest, parse_request
from services.seq_resolve import GeneResolver, SeqError
from services.seq_results import (
    purge_expired,
    read_status,
    result_dir,
    start_job,
    valid_id,
    write_status,
)
from services.seq_source import Source
from settings import settings

router = APIRouter()
logger = logging.getLogger(__name__)


class SequencesBody(BaseModel):
    """Exactly one of ``query``, ``loci`` or ``items`` selects the genes."""

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
        raise HTTPException(status_code=422, detail=exc.as_detail()) from exc


def _urls(request: Request, rid: str) -> dict[str, str]:
    base = (settings.PUBLIC_BASE_URL or str(request.base_url)).rstrip("/")
    return {"status_url": f"{base}/sequence-jobs/{rid}",
            "result_url": f"{base}/sequences/results/{rid}",
            "manifest_url": f"{base}/sequences/results/{rid}/manifest.tsv"}


def _open_source(req: SeqRequest) -> Source:
    try:
        return Source(req)
    except SeqError as exc:
        raise _http(exc) from exc
    except QueryError as exc:
        raise HTTPException(status_code=422, detail=exc.as_detail()) from exc


def _count(src: Source, req: SeqRequest):
    try:
        return src.count(req.max_genes)
    except QueryError as exc:
        raise HTTPException(status_code=422, detail=exc.as_detail()) from exc


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


def _job_payload(request: Request, rid: str, status: dict) -> dict:
    return {"job_id": rid, "result_id": rid, **_urls(request, rid),
            **{k: status.get(k) for k in ("state", "done", "total", "n_ok", "n_skipped",
                                          "n_sequences", "reasons", "error")}}


def _summary(result, data_version: str) -> dict:
    return {"n_sequences": result.n_sequences, "n_genes": result.n_genes,
            "n_skipped": result.n_skipped, "reasons": result.reasons,
            "data_version": data_version, "error": None}


@router.post("/sequences")
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
            return JSONResponse(_job_payload(request, rid, status), status_code=202,
                                headers={"Location": _urls(request, rid)["status_url"]})

        counted = _count(src, req)
        if not counted.complete:
            raise HTTPException(status_code=413, detail=error_detail(
                "too_many_genes", f"the request matches more than {req.max_genes} genes; narrow "
                "it or raise max_genes", n_genes_at_least=counted.n_genes, max_genes=req.max_genes))
        if counted.n_genes == 0:
            raise HTTPException(status_code=404, detail=error_detail(
                "no_genes_matched", "no gene matches the request", terms=src.terms()))
        data_version = src.data_version

        if counted.n_genes <= int(settings.SEQUENCES_SYNC_MAX_GENES):
            try:
                with sequences_slot():
                    write_status(rid, state="running", done=0, total=counted.n_genes)
                    t0 = time.monotonic()
                    try:
                        result = run_export(req, src.items(), out_dir=result_dir(rid),
                                            total=counted.n_genes)
                    except SeqError as exc:
                        shutil.rmtree(result_dir(rid), ignore_errors=True)
                        raise _http(exc) from exc
                    except QueryError as exc:
                        shutil.rmtree(result_dir(rid), ignore_errors=True)
                        raise HTTPException(status_code=422, detail=exc.as_detail()) from exc
                    status = write_status(rid, state="done", **_summary(result, data_version))
                    logger.info("POST /sequences sync genes=%d seqs=%d skipped=%d wall_s=%.2f",
                                result.n_genes, result.n_sequences, result.n_skipped,
                                time.monotonic() - t0)
            except BusyError:
                return _busy()
            return _fasta_response(request, rid, status)

        # Larger request: a job. It opens its own connections (this request's close below).
        total = counted.n_genes

        def work(out_dir, progress):
            job_src = Source(req)
            try:
                with sequences_slot(blocking=True):
                    result = run_export(req, job_src.items(), out_dir=out_dir, total=total,
                                        progress=progress)
                return _summary(result, job_src.data_version)
            finally:
                job_src.close()

        start_job(rid, work)
        status = read_status(rid) or {"state": "queued"}
        return JSONResponse(_job_payload(request, rid, status), status_code=202,
                            headers={"Location": _urls(request, rid)["status_url"]})
    finally:
        src.close()


@router.post("/sequences/estimate")
def post_sequences_estimate(body: SequencesBody, request: Request):
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
            "mode": "sync" if counted.n_genes <= int(settings.SEQUENCES_SYNC_MAX_GENES) else "job",
            "sync_max_genes": int(settings.SEQUENCES_SYNC_MAX_GENES),
            "cached": bool(status and status["state"] == "done"),
            "terms": src.terms(),
        }
    finally:
        src.close()


@router.post("/sequences/plan")
def post_sequences_plan(body: SequencesBody):
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
        raise HTTPException(status_code=422, detail=exc.as_detail()) from exc


@router.get("/sequence-jobs/{result_id}")
def get_sequence_job(result_id: str, request: Request):
    status = read_status(result_id) if valid_id(result_id) else None
    if status is None:
        raise HTTPException(status_code=404, detail=error_detail(
            "job_not_found", "no such job (results expire after a while)"))
    return JSONResponse(_job_payload(request, result_id, status),
                        headers={"Cache-Control": "no-store"})


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
    return status, result_dir(result_id)


@router.get("/sequences/results/{result_id}")
def get_sequence_result(result_id: str, request: Request):
    status, _ = _ready(result_id)
    return _fasta_response(request, result_id, status)


@router.get("/sequences/results/{result_id}/manifest.tsv")
def get_sequence_manifest(result_id: str):
    _status, d = _ready(result_id)
    return FileResponse(d / "manifest.tsv", media_type="text/tab-separated-values",
                        filename=f"annotrieve-{result_id[:8]}-manifest.tsv",
                        headers={"Cache-Control": "no-store"})
