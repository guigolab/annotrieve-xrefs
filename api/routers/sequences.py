"""Bulk transcript FASTA export (``POST /sequences``)."""
from __future__ import annotations

from typing import Literal

from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field

from helpers.concurrency import BusyError, sequences_slot
from helpers.deps import corpus_conn
from helpers.errors import error_detail, http_error
from helpers.pagination import (
    MAX_SEQUENCES_BUFFER_BP,
    MAX_SEQUENCES_ID_CHARS,
    MAX_SEQUENCES_POST,
)
from helpers.session import GeneCorpusUnavailable
from services.sequences import SequencesError, build_request, export_sequences

router = APIRouter()


class SequenceItemBody(BaseModel):
    annotation_id: str = Field(..., min_length=1, max_length=MAX_SEQUENCES_ID_CHARS)
    local_id: int = Field(..., ge=0)
    transcript_id: str = Field(..., min_length=1, max_length=MAX_SEQUENCES_ID_CHARS)


class SequencesBody(BaseModel):
    """JSON body for ``POST /sequences`` (bulk POST only)."""

    action: Literal["cds", "transcript", "protein"] = Field(
        ...,
        description=(
            "gffread extraction: cds (-x), transcript/spliced exons+UTRs (-w), "
            "protein (-y)"
        ),
    )
    items: list[SequenceItemBody] = Field(
        ...,
        min_length=1,
        description=(
            f"1..{MAX_SEQUENCES_POST} gene/transcript references "
            f"(cap enforced as too_many_items). "
            "transcript_id is required (no server-side longest-transcript default)."
        ),
    )
    buffer_bp: int = Field(
        0,
        ge=0,
        le=MAX_SEQUENCES_BUFFER_BP,
        description=(
            "Extra flank (bp) on each side of each transcript; passed to gffread "
            f"as --w-add. Only valid with action=transcript (0..{MAX_SEQUENCES_BUFFER_BP})."
        ),
    )
    include_exon_coords: bool = Field(
        False,
        description="Pass gffread -W (include exon coordinates on the defline).",
    )


def _busy_response() -> JSONResponse:
    return JSONResponse(
        status_code=503,
        content={
            "detail": error_detail(
                "sequences_busy",
                "too many concurrent sequences requests; retry shortly",
            )
        },
        headers={"Retry-After": "1"},
    )


def _http_for_sequences_error(exc: SequencesError) -> HTTPException:
    headers = (
        {"Retry-After": str(exc.retry_after)}
        if exc.retry_after is not None
        else None
    )
    return HTTPException(
        status_code=exc.status,
        detail=error_detail(exc.code, str(exc), **exc.extra),
        headers=headers,
    )


def _cleanup_temp(temp_dir) -> None:
    try:
        temp_dir.cleanup()
    except Exception:
        pass


@router.post("/sequences")
def post_sequences(body: SequencesBody, background_tasks: BackgroundTasks):
    """
    Export a multi-FASTA for the requested transcripts.

    Fail-first: any unresolved seqid / missing transcript / NCBI error aborts
    with no partial body. Response is ``text/x-fasta``.
    """
    try:
        request = build_request(
            action=body.action,
            items=body.items,
            buffer_bp=body.buffer_bp,
            include_exon_coords=body.include_exon_coords,
        )
    except SequencesError as exc:
        raise _http_for_sequences_error(exc) from exc

    try:
        with sequences_slot():
            corpus = corpus_conn()
            try:
                result = export_sequences(request, corpus=corpus)
            except GeneCorpusUnavailable as exc:
                raise http_error(
                    503,
                    "gene_corpus_unavailable",
                    str(exc),
                ) from exc
            except SequencesError as exc:
                raise _http_for_sequences_error(exc) from exc
    except BusyError:
        return _busy_response()

    background_tasks.add_task(_cleanup_temp, result.temp_dir)

    return FileResponse(
        path=result.fasta_path,
        media_type="text/x-fasta",
        filename=result.filename,
        headers={
            "X-Sequence-Count": str(result.sequence_count),
            "Cache-Control": "no-store",
        },
    )
