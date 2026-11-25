"""FASTA downloads (``/downloads``) and transcript choice (``/transcripts``, ``/genes/{ref}/transcripts``).

``POST /downloads`` puts a request in the one queue and answers 202 with a link; a runner (one at a
time for the whole container, see ``seq_runner``) produces a gzipped FASTA and a report. The id of
a download is a hash of the normalised request and the data version: an equal request finds the
download that exists instead of making a second one.
"""
from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from typing import Any

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel, ConfigDict, Field

from helpers.client_ip import client_ip, owner_key
from helpers.concurrency import plan_slot
from helpers.cursor import QueryError
from helpers.errors import error_detail, http_for_query_error
from services.seq_export import count_sequences, run_export
from services.seq_request import SeqRequest, parse_request
from services.seq_resolve import GeneResolver, SeqError
from services.seq_results import (
    entries,
    Refused,
    admission_lock,
    admit,
    enqueue,
    fasta_path,
    mark_downloaded,
    purge_expired,
    quota_status,
    queue_position,
    read_status,
    remove,
    result_dir,
    valid_id,
)
from services.seq_source import Source
from settings import settings

router = APIRouter()

_ERR = {"description": "Error: `detail` is `{code, message, ...}`."}
_FASTA_HEADERS = {
    "X-Sequence-Count": {"description": "Sequences in the FASTA", "schema": {"type": "integer"}},
    "X-Skipped-Count": {"description": "Genes skipped (see the report)", "schema": {"type": "integer"}},
    "X-Gene-Count": {"description": "Genes processed", "schema": {"type": "integer"}},
    "X-Download-Id": {"description": "Identifier of this download", "schema": {"type": "string"}},
    "X-Data-Version": {"description": "Index version used", "schema": {"type": "string"}},
}
POST_RESPONSES = {
    200: {"description": "`dry_run=true`: the estimate. Otherwise the finished download already "
                         "stored for an equal request (same body as GET /downloads/{id})."},
    202: {"description": "Queued (or already queued/running for an equal request). "
                         "`Location` is the status URL; poll it with `Retry-After`."},
    404: {**_ERR, "description": "`no_genes_matched`, or a pinned/unknown locus (loci mode, on_error=fail)"},
    413: {**_ERR, "description": "`too_many_genes` or `too_many_sequences` (over the limit of one download)"},
    422: {**_ERR, "description": "Invalid request (`validation_error`, `exactly_one_source`, "
                                 "`unknown_rank_key`, `invalid_locus_ref`, ...)"},
    429: {**_ERR, "description": "`too_many_active_downloads` or `daily_quota_exceeded` for this "
                                 "client IP (see `Retry-After` and GET /downloads/quota)"},
    503: {**_ERR, "description": "`downloads_queue_full`, `downloads_capacity`, "
                                 "`annotations_unavailable`, `locus_ref_index_unavailable`, "
                                 "`plan_busy` (too many GFF-reading requests at once; see Retry-After)"},
}


logger = logging.getLogger(__name__)


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


class TranscriptPlan(BaseModel):
    """The transcripts a request would extract, with the metrics behind the choice."""

    data_version: str
    action: str
    n_genes: int
    n_skipped: int
    genes: list[PlanGene]
    loci: list[PlanLocus] = Field(description="the selection, ready to post to /downloads as `loci`")


class DownloadBody(BaseModel):
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


def _parse(body: DownloadBody, *, default_action: str | None = None) -> SeqRequest:
    data = body.model_dump(exclude_none=True)
    if default_action and "action" not in data:
        data["action"] = default_action
    try:
        return parse_request(data)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc


def _base(request: Request) -> str:
    return (settings.PUBLIC_BASE_URL or str(request.base_url)).rstrip("/")


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


def _payload(request: Request, rid: str, st: dict) -> dict:
    """The status of a download, as every route shows it."""
    base = f"{_base(request)}/downloads/{rid}"
    state = st.get("state")
    links = {"self": base}
    if state == "done":
        links.update(fasta=f"{base}/fasta", report=f"{base}/report")
    out: dict[str, Any] = {
        "id": rid, "output": "fasta", "state": state,
        "created_at": st.get("created"),
        "queue_position": queue_position(rid) if state == "queued" else None,
        "progress": {"genes_done": st.get("done", 0), "genes_total": st.get("total")},
        "n_sequences": st.get("n_sequences"), "n_genes": st.get("n_genes"),
        "n_skipped": st.get("n_skipped"), "reasons": st.get("reasons"),
        "size_bytes": st.get("size_bytes"), "data_version": st.get("data_version"),
        "error": st.get("error"), "links": links,
    }
    if state == "done":
        keep = (int(settings.DOWNLOAD_KEEP_AFTER_FETCH_S) if st.get("downloaded_at")
                else int(settings.DOWNLOAD_KEEP_UNFETCHED_S))
        out["expires_at"] = (st.get("downloaded_at") or st.get("updated", 0)) + keep
    return out


def _reply(request: Request, rid: str, st: dict, status_code: int) -> JSONResponse:
    headers = {"Location": f"{_base(request)}/downloads/{rid}"}
    if st.get("state") in ("queued", "running"):
        headers["Retry-After"] = "5"
    return JSONResponse(_payload(request, rid, st), status_code=status_code, headers=headers)


def _stored_file(rid: str) -> bool:
    return fasta_path(rid).is_file()


@router.get("/downloads/quota", summary="Downloads left for this client today",
            responses={200: {"description": "Rolling 24-hour quota of the calling IP."}})
def get_download_quota(request: Request):
    """New downloads used and left in the last 24 hours (an equal request is free), and the
    downloads of this client that are waiting or running now."""
    owner = owner_key(client_ip(request))
    q = quota_status(owner)
    active = sum(1 for _d, st in entries() if st and st.get("owner") == owner
                 and st.get("state") in ("queued", "running"))
    return {**q, "active": active, "max_active": int(settings.SEQUENCES_IP_MAX_ACTIVE)}


@router.post("/downloads", responses=POST_RESPONSES, summary="Ask for a FASTA download",
             response_model=None)
def create_download(body: DownloadBody, request: Request,
                    dry_run: bool = Query(False, description="only count and estimate; nothing is "
                                                             "queued and no quota is used")):
    """Queue a CDS, transcript or protein FASTA for the genes of a query or an explicit list.

    All downloads wait in one queue (one runs at a time, because of the NCBI request budget). The
    answer is `202` with the status; poll `GET /downloads/{id}` and fetch `/fasta` when `state`
    is `done`. An equal request returns the download that exists.
    """
    req = _parse(body, default_action="cds" if dry_run else None)
    src = _open_source(req)
    try:
        rid = req.result_id(src.data_version)
        purge_expired()
        status = read_status(rid)
        if dry_run:
            counted = _count(src, req)
            return {
                "data_version": src.data_version, "id": rid,
                "n_genes": counted.n_genes if counted.complete else None,
                "n_genes_at_least": counted.n_genes, "complete": counted.complete,
                "est_bp_upper_bound": counted.bp_upper_bound, "max_genes": req.max_genes,
                "over_limit": not counted.complete,
                "max_sequences": int(settings.SEQUENCES_MAX_SEQUENCES),
                "cached": bool(status and status["state"] == "done" and _stored_file(rid)),
                "terms": src.terms(),
            }
        if status and status["state"] == "done" and _stored_file(rid):
            return _reply(request, rid, status, 200)
        if status and status["state"] in ("queued", "running"):
            return _reply(request, rid, status, 202)

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
            # refuse before queueing.
            limit = int(settings.SEQUENCES_MAX_SEQUENCES)
            with plan_slot():
                n_seq = count_sequences(req, src.items(), limit)
            if n_seq > limit:
                raise HTTPException(status_code=413, detail=error_detail(
                    "too_many_sequences", f"the request gives more than {limit} sequences; narrow it "
                    "or choose one transcript per gene", n_sequences_at_least=n_seq,
                    max_sequences=limit))

        owner = owner_key(client_ip(request))
        with admission_lock():
            # Look again now that no other request can change the state: an equal request that
            # arrived while this one was counting must not be queued twice or charged twice.
            status = read_status(rid)
            if status and status["state"] == "done" and _stored_file(rid):
                return _reply(request, rid, status, 200)
            if status and status["state"] in ("queued", "running"):
                return _reply(request, rid, status, 202)
            try:
                admit(owner)
            except Refused as exc:
                raise HTTPException(
                    status_code=exc.status, detail=error_detail(exc.code, str(exc), **exc.extra),
                    headers={"Retry-After": str(exc.retry_after)}) from exc
            stored = enqueue(rid, body.model_dump(exclude_none=True) | {"action": req.action.value},
                             owner=owner, total=counted.n_genes, data_version=src.data_version)
        return _reply(request, rid, stored, 202)
    finally:
        src.close()


def _build_plan(req: SeqRequest, src: Source, *, cap: int) -> dict:
    counted = _count_cap(src, cap)
    if not counted.complete:
        raise HTTPException(status_code=413, detail=error_detail(
            "plan_too_large", f"this covers at most {cap} genes; narrow the request",
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


@router.post("/transcripts", response_model=TranscriptPlan,
             summary="Which transcripts a request would pick",
             responses={413: _ERR, 422: _ERR, 503: _ERR})
def list_transcripts(body: DownloadBody):
    """The candidate transcripts of every gene of the request, with the metrics behind the choice
    and the ones the rule selects. Same body as `POST /downloads`.

    Reads the GFF slices (no network, nothing stored). A client can change the selection and post
    the resulting ``loci`` to ``/downloads``.
    """
    req = _parse(body)
    src = _open_source(req)
    try:
        with plan_slot():
            return _build_plan(req, src, cap=min(req.max_genes, int(settings.SEQUENCES_PLAN_MAX_GENES)))
    finally:
        src.close()


@router.get("/genes/{locus_ref}/transcripts", response_model=PlanGene,
            summary="Transcripts of one gene", responses={404: _ERR, 422: _ERR, 503: _ERR})
def get_gene_transcripts(
    locus_ref: str,
    action: str = Query("cds", description="cds | transcript | protein: decides which transcripts "
                                           "are eligible (`cds`/`protein` need a CDS)"),
    transcript: str | None = Query(None, description="selection preset (default `longest`)"),
):
    """Every transcript of one gene (`annotation_id:local_id`) with its metrics, and which one
    the preset would select."""
    data: dict[str, Any] = {"action": action, "loci": [locus_ref]}
    if transcript:
        data["transcript"] = transcript
    try:
        req = parse_request(data)
    except QueryError as exc:
        raise http_for_query_error(exc) from exc
    src = _open_source(req)
    try:
        with plan_slot():
            gene = _build_plan(req, src, cap=1)["genes"][0]
    finally:
        src.close()
    if gene["status"] == "skipped" and gene["reason"] in ("annotation_not_found", "locus_not_found"):
        raise HTTPException(status_code=404, detail=error_detail(gene["reason"], gene["message"]))
    return gene


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


def _known(rid: str) -> dict:
    st = read_status(rid) if valid_id(rid) else None
    if st is None:
        raise HTTPException(status_code=404, detail=error_detail(
            "download_not_found", "no such download (downloads are removed after a while)"))
    return st


@router.get("/downloads/{download_id}", summary="Status of a download",
            responses={404: _ERR})
def get_download(download_id: str, request: Request):
    """`state` is `queued`, `running`, `done` or `failed`. While waiting, `queue_position` tells how
    many downloads are ahead; while running, `progress` counts genes. `Retry-After` says when to
    ask again."""
    st = _known(download_id)
    headers = {"Retry-After": "5"} if st.get("state") in ("queued", "running") else None
    return JSONResponse(_payload(request, download_id, st), headers=headers)


def _ready(rid: str) -> dict:
    st = _known(rid)
    if st["state"] == "failed":
        err = st.get("error") or {}
        raise HTTPException(status_code=int(err.get("status", 500)), detail=error_detail(
            err.get("code", "download_failed"), err.get("message", "the download failed")))
    if st["state"] != "done":
        raise HTTPException(status_code=409, detail=error_detail(
            "download_not_ready", "the download has not finished", state=st["state"],
            queue_position=queue_position(rid), done=st.get("done"), total=st.get("total")),
            headers={"Retry-After": "5"})
    if not _stored_file(rid):
        raise HTTPException(status_code=404, detail=error_detail(
            "download_not_found", "this download has been removed; ask for it again"))
    return st


def _stream(fh, rid: str) -> Iterator[bytes]:
    """The open file in blocks; once all of it has been sent the short keep-time starts."""
    with fh:
        while block := fh.read(1024 * 1024):
            yield block
    mark_downloaded(rid)


@router.get("/downloads/{download_id}/fasta", summary="The FASTA (gzip)",
            responses={200: {"description": "gzip-compressed FASTA (`application/gzip`). "
                                            "Removed 10 minutes after the first complete download.",
                             "headers": _FASTA_HEADERS,
                             "content": {"application/gzip": {"schema": {"type": "string",
                                                                          "format": "binary"}}}},
                       404: _ERR, 409: {**_ERR, "description": "`download_not_ready`"}})
def get_download_fasta(download_id: str):
    st = _ready(download_id)
    try:
        fh = open(fasta_path(download_id), "rb")        # opened now: a purge cannot unlink it half-way
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=error_detail(
            "download_not_found", "this download has been removed; ask for it again")) from None
    size = os.fstat(fh.fileno()).st_size
    return StreamingResponse(
        _stream(fh, download_id), media_type="application/gzip",
        headers={"Content-Disposition": f'attachment; filename="annotrieve-{download_id[:8]}.fa.gz"',
                 "Content-Length": str(size), "Cache-Control": "no-store",
                 "X-Sequence-Count": str(st.get("n_sequences", 0)),
                 "X-Skipped-Count": str(st.get("n_skipped", 0)),
                 "X-Gene-Count": str(st.get("n_genes", 0)),
                 "X-Download-Id": download_id, "X-Data-Version": str(st.get("data_version", ""))})


@router.get("/downloads/{download_id}/report", summary="Per-gene outcome table (tsv)",
            responses={200: {"description": "tsv", "content": {"text/tab-separated-values": {
                "schema": {"type": "string"}}}}, 404: _ERR, 409: _ERR})
def get_download_report(download_id: str):
    _ready(download_id)
    manifest = result_dir(download_id) / "manifest.tsv"
    if not manifest.is_file():
        raise HTTPException(status_code=404, detail=error_detail(
            "download_not_found", "this download has been removed; ask for it again"))
    return FileResponse(manifest,
                        media_type="text/tab-separated-values",
                        filename=f"annotrieve-{download_id[:8]}-report.tsv",
                        headers={"Cache-Control": "no-store"})


@router.delete("/downloads/{download_id}", status_code=204, summary="Cancel or delete a download",
               responses={404: _ERR})
def delete_download(download_id: str, request: Request):
    """Takes it out of the queue, stops it if running (at the next batch), or deletes the files.
    Only the client that asked for it may delete it."""
    st = _known(download_id)
    if st.get("owner") != owner_key(client_ip(request)):
        raise HTTPException(status_code=403, detail=error_detail(
            "not_your_download", "only the client that asked for this download can delete it"))
    if st.get("state") == "running":
        (result_dir(download_id) / "cancel").touch()
    else:
        remove(download_id)
    from fastapi import Response
    return Response(status_code=204)
