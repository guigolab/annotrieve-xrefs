"""Bulk transcript FASTA export (POST /sequences).

Items are grouped by ``(annotation_id, local_id)`` so each gene costs one
tabix query and one gffread run. NCBI windows on the same accession that
overlap or sit within ``SEQUENCES_COALESCE_GAP_BP`` are merged into a single
efetch (fewer round-trips), then sliced per gene in memory.
"""
from __future__ import annotations

import logging
import shutil
import sqlite3
import threading
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Iterable

import pysam

from helpers.annotation_files import (
    AnnotationFiles,
    AnnotationFilesError,
    resolve_annotation_files,
    resolve_seqid_accession,
)
from helpers.gff_slice import (
    TranscriptModel,
    dedupe_gff_lines,
    load_features_from_lines,
    resolve_transcripts,
    shift_gff_lines,
    write_fasta,
    write_gff,
)
from helpers.gffread_runner import (
    GffreadError,
    SequenceAction,
    count_fasta_records,
    run_gffread,
    split_fasta_to_record_files,
)
from helpers.ncbi_eutils import NcbiError, get_ncbi_client
from helpers.pagination import (
    MAX_SEQUENCES_BUFFER_BP,
    MAX_SEQUENCES_ID_CHARS,
    MAX_SEQUENCES_POST,
)
from helpers.windows import GenomicWindow, coalesce_windows, slice_sequence
from services.annotations import ShardUnavailable, get_gene, open_shard
from settings import settings

logger = logging.getLogger(__name__)

_TABIX_LOCAL = threading.local()
_TABIX_MAX_HANDLES = 32
_TABIX_REGISTRY_LOCK = threading.Lock()
# thread ident -> OrderedDict of path -> TabixFile (for process-wide close)
_TABIX_REGISTRY: dict[int, OrderedDict[str, pysam.TabixFile]] = {}


class SequencesError(Exception):
    """Fail-first sequences pipeline error."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        status: int = 422,
        retry_after: str | None = None,
        **extra: Any,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.status = status
        self.retry_after = retry_after
        self.extra = extra


@dataclass(frozen=True)
class SequenceItem:
    annotation_id: str
    local_id: int
    transcript_id: str
    index: int  # position in the deduped request list


@dataclass
class GeneGroup:
    annotation_id: str
    local_id: int
    items: list[SequenceItem] = field(default_factory=list)
    # Filled during resolve:
    gene: dict | None = None
    annotation: dict | None = None
    files: AnnotationFiles | None = None
    accession: str | None = None
    transcripts: dict[str, TranscriptModel] = field(default_factory=dict)
    window: GenomicWindow | None = None
    group_key: str = ""


@dataclass(frozen=True)
class SequencesRequest:
    action: SequenceAction
    items: list[SequenceItem]
    buffer_bp: int = 0
    include_exon_coords: bool = False


@dataclass
class SequencesResult:
    fasta_path: Path
    sequence_count: int
    temp_dir: TemporaryDirectory
    filename: str


def _validate_id(value: str, *, field: str) -> str:
    text = (value or "").strip()
    if not text:
        raise SequencesError(
            "validation_error",
            f"{field} is required",
            status=422,
            field=field,
        )
    if len(text) > MAX_SEQUENCES_ID_CHARS:
        raise SequencesError(
            "validation_error",
            f"{field} exceeds {MAX_SEQUENCES_ID_CHARS} characters",
            status=422,
            field=field,
        )
    if any(ord(ch) < 32 for ch in text):
        raise SequencesError(
            "validation_error",
            f"{field} contains control characters",
            status=422,
            field=field,
        )
    return text


def parse_sequence_items(raw_items: Iterable[dict | Any]) -> list[SequenceItem]:
    """Validate / dedupe items preserving first-seen order."""
    items_list = list(raw_items or [])
    if not items_list:
        raise SequencesError(
            "items_required",
            "items must contain at least one entry",
            status=422,
        )
    if len(items_list) > MAX_SEQUENCES_POST:
        raise SequencesError(
            "too_many_items",
            f"at most {MAX_SEQUENCES_POST} items allowed",
            status=422,
            max_items=MAX_SEQUENCES_POST,
            got=len(items_list),
        )

    seen: set[tuple[str, int, str]] = set()
    out: list[SequenceItem] = []
    for i, raw in enumerate(items_list):
        if hasattr(raw, "model_dump"):
            data = raw.model_dump()
        elif isinstance(raw, dict):
            data = raw
        else:
            data = {
                "annotation_id": getattr(raw, "annotation_id", None),
                "local_id": getattr(raw, "local_id", None),
                "transcript_id": getattr(raw, "transcript_id", None),
            }
        annotation_id = _validate_id(
            str(data.get("annotation_id") or ""),
            field="annotation_id",
        )
        transcript_id = _validate_id(
            str(data.get("transcript_id") or ""),
            field="transcript_id",
        )
        try:
            local_id = int(data.get("local_id"))
        except (TypeError, ValueError) as exc:
            raise SequencesError(
                "validation_error",
                "local_id must be an integer",
                status=422,
                index=i,
            ) from exc
        if local_id < 0:
            raise SequencesError(
                "validation_error",
                "local_id must be >= 0",
                status=422,
                index=i,
            )
        key = (annotation_id, local_id, transcript_id)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            SequenceItem(
                annotation_id=annotation_id,
                local_id=local_id,
                transcript_id=transcript_id,
                index=len(out),
            )
        )
    if not out:
        raise SequencesError(
            "items_required",
            "items must contain at least one entry",
            status=422,
        )
    return out


def build_request(
    *,
    action: str | SequenceAction,
    items: Iterable[dict | Any],
    buffer_bp: int = 0,
    include_exon_coords: bool = False,
) -> SequencesRequest:
    try:
        action_enum = (
            action if isinstance(action, SequenceAction) else SequenceAction(action)
        )
    except ValueError as exc:
        raise SequencesError(
            "validation_error",
            f"unsupported action {action!r}",
            status=422,
        ) from exc

    try:
        buf = int(buffer_bp)
    except (TypeError, ValueError) as exc:
        raise SequencesError(
            "validation_error",
            "buffer_bp must be an integer",
            status=422,
        ) from exc
    if buf < 0 or buf > MAX_SEQUENCES_BUFFER_BP:
        raise SequencesError(
            "validation_error",
            f"buffer_bp must be between 0 and {MAX_SEQUENCES_BUFFER_BP}",
            status=422,
            max_buffer_bp=MAX_SEQUENCES_BUFFER_BP,
            buffer_bp=buf,
        )
    if buf > 0 and action_enum != SequenceAction.transcript:
        raise SequencesError(
            "buffer_not_supported",
            "buffer_bp is only supported with action=transcript",
            status=422,
            action=action_enum.value,
            buffer_bp=buf,
        )

    return SequencesRequest(
        action=action_enum,
        items=parse_sequence_items(items),
        buffer_bp=buf,
        include_exon_coords=bool(include_exon_coords),
    )


def _deadline_check(deadline: float) -> None:
    if time.monotonic() > deadline:
        raise SequencesError(
            "sequences_timeout",
            "sequences request exceeded deadline",
            status=504,
        )


def _annotation_row(
    corpus: sqlite3.Connection,
    annotation_id: str,
) -> dict:
    row = corpus.execute(
        """
        SELECT annotation_id, taxid, assembly_accession, organism_name,
               source_database, profile_id, gff_path
        FROM annotation
        WHERE annotation_id = ?
        """,
        (annotation_id,),
    ).fetchone()
    if row is None:
        raise SequencesError(
            "annotation_not_found",
            f"annotation {annotation_id} not found",
            status=404,
            annotation_id=annotation_id,
        )
    return {
        "annotation_id": str(row["annotation_id"]),
        "taxid": int(row["taxid"]),
        "assembly_accession": str(row["assembly_accession"]),
        "organism_name": row["organism_name"],
        "source_database": str(row["source_database"] or ""),
        "profile_id": str(row["profile_id"] or ""),
        "gff_path": row["gff_path"] or "",
    }


def _open_tabix(gff_path: Path, index_path: Path | None = None):
    """Thread-local LRU of TabixFile handles (pysam is not thread-safe).

    Annotrieve indexes with ``tabix --csi`` (``.csi``), not ``.tbi``. Pass the
    resolved index path so pysam does not default to looking for ``.tbi``.
    """
    cache: OrderedDict[str, pysam.TabixFile] | None = getattr(
        _TABIX_LOCAL, "handles", None
    )
    if cache is None:
        cache = OrderedDict()
        _TABIX_LOCAL.handles = cache
        with _TABIX_REGISTRY_LOCK:
            _TABIX_REGISTRY[threading.get_ident()] = cache
    index = str(index_path) if index_path is not None else None
    key = f"{gff_path}\0{index or ''}"
    handle = cache.get(key)
    if handle is not None:
        cache.move_to_end(key)
        return handle
    handle = (
        pysam.TabixFile(str(gff_path), index=index)
        if index
        else pysam.TabixFile(str(gff_path))
    )
    cache[key] = handle
    while len(cache) > _TABIX_MAX_HANDLES:
        _old_key, old = cache.popitem(last=False)
        try:
            old.close()
        except Exception:
            pass
    return handle


def close_tabix_handles() -> None:
    """Close all thread-local tabix handles (tests / shutdown)."""
    with _TABIX_REGISTRY_LOCK:
        caches = list(_TABIX_REGISTRY.values())
        _TABIX_REGISTRY.clear()
    for cache in caches:
        while cache:
            _key, handle = cache.popitem(last=False)
            try:
                handle.close()
            except Exception:
                pass
    if getattr(_TABIX_LOCAL, "handles", None) is not None:
        _TABIX_LOCAL.handles = None


def _fetch_gene_lines(
    files: AnnotationFiles,
    *,
    seqid: str,
    gene_start: int,
    gene_end: int,
) -> list[str]:
    # pysam regions are 0-based half-open.
    start0 = max(0, int(gene_start) - 1)
    end0 = max(start0 + 1, int(gene_end))
    try:
        tbx = _open_tabix(files.gff_path, files.csi_path)
        return list(tbx.fetch(seqid, start0, end0))
    except ValueError as exc:
        # Unknown contig etc.
        raise SequencesError(
            "seqid_not_resolved",
            f"tabix fetch failed for seqid {seqid!r}: {exc}",
            status=422,
            seqid=seqid,
        ) from exc
    except OSError as exc:
        raise SequencesError(
            "annotation_file_not_found",
            f"failed to read GFF: {exc}",
            status=404,
            gff_path=str(files.gff_path),
        ) from exc
    except Exception as exc:
        # pysam.SamtoolsError and other htslib failures are bare Exception.
        raise SequencesError(
            "annotation_file_not_found",
            f"failed to read GFF: {exc}",
            status=404,
            gff_path=str(files.gff_path),
        ) from exc


def _group_items(items: list[SequenceItem]) -> list[GeneGroup]:
    groups: dict[tuple[str, int], GeneGroup] = {}
    order: list[tuple[str, int]] = []
    for item in items:
        key = (item.annotation_id, item.local_id)
        if key not in groups:
            groups[key] = GeneGroup(
                annotation_id=item.annotation_id,
                local_id=item.local_id,
                group_key=f"{item.annotation_id}:{item.local_id}",
            )
            order.append(key)
        groups[key].items.append(item)
    return [groups[k] for k in order]


def _item_extra(item: SequenceItem, **more: Any) -> dict[str, Any]:
    out = {
        "index": item.index,
        "annotation_id": item.annotation_id,
        "local_id": item.local_id,
        "transcript_id": item.transcript_id,
    }
    out.update(more)
    return out


_RESERVED_EXTRA_KEYS = frozenset({"status", "retry_after", "code", "message"})


def _merge_extra(*parts: dict[str, Any]) -> dict[str, Any]:
    """Merge error extras; drop keys that collide with SequencesError kwargs."""
    out: dict[str, Any] = {}
    for part in parts:
        if not part:
            continue
        for key, value in part.items():
            if key in _RESERVED_EXTRA_KEYS:
                continue
            out[key] = value
    return out


def _resolve_group(
    group: GeneGroup,
    *,
    corpus: sqlite3.Connection,
    shard_cache: dict[str, sqlite3.Connection],
    ann_cache: dict[str, dict],
    files_cache: dict[str, AnnotationFiles],
    accession_cache: dict[tuple[str, str], str],
    action: SequenceAction,
    buffer_bp: int,
) -> None:
    first = group.items[0]
    try:
        if group.annotation_id not in ann_cache:
            ann_cache[group.annotation_id] = _annotation_row(
                corpus, group.annotation_id
            )
        annotation = ann_cache[group.annotation_id]
        group.annotation = annotation

        if group.annotation_id not in shard_cache:
            try:
                shard_cache[group.annotation_id] = open_shard(group.annotation_id)
            except ShardUnavailable as exc:
                raise SequencesError(
                    "shard_not_found",
                    str(exc),
                    status=404,
                    annotation_id=group.annotation_id,
                ) from exc
        shard = shard_cache[group.annotation_id]
        gene = get_gene(shard, group.local_id)
        if gene is None:
            raise SequencesError(
                "gene_not_found",
                f"gene local_id {group.local_id} not found",
                status=404,
                **_item_extra(first),
            )
        group.gene = gene

        if group.annotation_id not in files_cache:
            try:
                files_cache[group.annotation_id] = resolve_annotation_files(
                    taxid=int(annotation["taxid"]),
                    assembly_accession=str(annotation["assembly_accession"]),
                    gff_path=str(annotation["gff_path"] or ""),
                    source_database=str(annotation["source_database"] or ""),
                )
            except AnnotationFilesError as exc:
                status = 503 if exc.code == "annotations_unavailable" else 404
                if exc.code == "seqid_not_resolved":
                    status = 422
                raise SequencesError(
                    exc.code,
                    str(exc),
                    status=status,
                    **_merge_extra(_item_extra(first), exc.extra),
                ) from exc
        files = files_cache[group.annotation_id]
        group.files = files

        seqid = str(gene["seqid"])
        acc_key = (group.annotation_id, seqid)
        if acc_key not in accession_cache:
            try:
                accession_cache[acc_key] = resolve_seqid_accession(files, seqid)
            except AnnotationFilesError as exc:
                raise SequencesError(
                    exc.code,
                    str(exc),
                    status=422 if exc.code == "seqid_not_resolved" else 404,
                    **_merge_extra(
                        _item_extra(first, seqid=seqid),
                        exc.extra,
                    ),
                ) from exc
        group.accession = accession_cache[acc_key]

        lines = _fetch_gene_lines(
            files,
            seqid=seqid,
            gene_start=int(gene["start"]),
            gene_end=int(gene["end"]),
        )
        features = load_features_from_lines(lines)
        wanted = [it.transcript_id for it in group.items]
        try:
            models = resolve_transcripts(
                features,
                source_gene_id=str(gene["source_gene_id"]),
                transcript_ids=wanted,
            )
        except KeyError as exc:
            missing = str(exc.args[0]) if exc.args else "unknown"
            bad = next(
                (it for it in group.items if it.transcript_id == missing),
                first,
            )
            raise SequencesError(
                "transcript_not_found",
                f"transcript {missing!r} not found in gene window",
                status=404,
                **_item_extra(bad, seqid=seqid),
            ) from exc

        gene_start = int(gene["start"])
        gene_end = int(gene["end"])
        for item in group.items:
            model = models[item.transcript_id]
            strand = (model.strand or "").strip()
            if strand not in {"+", "-"}:
                raise SequencesError(
                    "transcript_strand_unknown",
                    f"transcript {item.transcript_id!r} has unknown strand",
                    status=422,
                    **_item_extra(item, seqid=seqid, strand=strand or "."),
                )
            if not model.has_structure:
                raise SequencesError(
                    "transcript_incomplete",
                    (
                        f"transcript {item.transcript_id!r} has no structural "
                        "children in the gene window (exons/CDS may fall "
                        "outside shard gene bounds)"
                    ),
                    status=422,
                    **_item_extra(item, seqid=seqid),
                )
            if model.start < gene_start or model.end > gene_end:
                raise SequencesError(
                    "transcript_outside_gene",
                    (
                        f"transcript {item.transcript_id!r} span "
                        f"{model.start}-{model.end} falls outside shard gene "
                        f"bounds {gene_start}-{gene_end}"
                    ),
                    status=422,
                    **_merge_extra(
                        _item_extra(item, seqid=seqid),
                        {
                            "gene_start": gene_start,
                            "gene_end": gene_end,
                            "tx_start": model.start,
                            "tx_end": model.end,
                        },
                    ),
                )
            if action in {SequenceAction.cds, SequenceAction.protein} and not model.has_cds:
                raise SequencesError(
                    "no_cds",
                    f"transcript {item.transcript_id!r} has no CDS",
                    status=422,
                    **_item_extra(item, seqid=seqid),
                )

        group.transcripts = models
        tx_start = min(m.start for m in models.values())
        tx_end = max(m.end for m in models.values())
        fetch_start = max(1, tx_start - buffer_bp)
        fetch_stop = tx_end + buffer_bp
        window_len = fetch_stop - fetch_start + 1
        max_window = int(settings.SEQUENCES_MAX_WINDOW_BP)
        if window_len > max_window:
            raise SequencesError(
                "window_too_large",
                f"fetch window {window_len} bp exceeds limit {max_window}",
                status=422,
                **_merge_extra(
                    _item_extra(first, seqid=seqid),
                    {
                        "window_bp": window_len,
                        "max_window_bp": max_window,
                    },
                ),
            )
        group.window = GenomicWindow(
            accession=group.accession,
            start=fetch_start,
            stop=fetch_stop,
        )
    except SequencesError:
        raise


def _species_label(annotation: dict) -> str:
    name = (annotation.get("organism_name") or "").strip()
    if name:
        return name.replace(" ", "_")
    return str(annotation.get("taxid") or "unknown")


def _defline(
    *,
    transcript_id: str,
    gene_name: str,
    species: str,
    annotation_id: str,
    seqid: str,
    tx_start: int,
    tx_end: int,
    strand: str,
    buffer_bp: int,
) -> str:
    short = annotation_id[:8]
    gene = gene_name or "gene"
    return (
        f"{transcript_id}|{gene}|{species}|{short} "
        f"{seqid}:{tx_start}-{tx_end}({strand}) buffer={buffer_bp}"
    )


def _group_error_extra(group: GeneGroup) -> dict[str, Any]:
    first = group.items[0]
    return {
        **_item_extra(first),
        "transcript_ids": [it.transcript_id for it in group.items],
    }


def _run_group_gffread(
    group: GeneGroup,
    *,
    sequence: str,
    window: GenomicWindow,
    action: SequenceAction,
    buffer_bp: int,
    include_exon_coords: bool,
    work_dir: Path,
    gid: str,
    deadline: float,
) -> dict[str, Path]:
    """
    Run gffread for one gene group.

    Returns ``{transcript_id: single-record FASTA path}``.
    """
    assert group.gene is not None
    assert group.annotation is not None
    assert group.files is not None

    # Collect GFF lines for requested transcripts, shifted into window coords.
    raw_lines: list[str] = []
    header_for_id: dict[str, str] = {}
    gene_name = str(
        group.gene.get("primary_name")
        or group.gene.get("source_gene_id")
        or "gene"
    )
    species = _species_label(group.annotation)
    seqid = str(group.gene["seqid"])

    for item in group.items:
        model = group.transcripts[item.transcript_id]
        raw_lines.extend(model.lines)
        header_for_id[item.transcript_id] = _defline(
            transcript_id=item.transcript_id,
            gene_name=gene_name,
            species=species,
            annotation_id=group.annotation_id,
            seqid=seqid,
            tx_start=model.start,
            tx_end=model.end,
            strand=model.strand,
            buffer_bp=buffer_bp,
        )

    raw_lines = dedupe_gff_lines(raw_lines)

    try:
        shifted = shift_gff_lines(raw_lines, new_seqid=gid, origin=window.start)
    except ValueError as exc:
        raise SequencesError(
            "gffread_failed",
            f"failed to build gffread input: {exc}",
            status=500,
            **_group_error_extra(group),
        ) from exc

    gff_path = work_dir / f"{gid}.gff3"
    fa_path = work_dir / f"{gid}.fa"
    out_raw = work_dir / f"{gid}.{action.value}.raw.fa"
    records_dir = work_dir / f"{gid}.records"
    write_gff(gff_path, shifted)
    write_fasta(fa_path, seqid=gid, sequence=sequence)

    # Clamp --w-add when the FASTA is shorter than the requested buffer flank
    # (chromosome start/end). gffread adds flank relative to the transcript
    # inside the provided FASTA.
    effective_buffer = buffer_bp
    if buffer_bp and action == SequenceAction.transcript:
        tx_start = min(m.start for m in group.transcripts.values())
        tx_end = max(m.end for m in group.transcripts.values())
        left_flank = tx_start - window.start
        right_flank = window.stop - tx_end
        seq_right = window.start + len(sequence) - 1
        if seq_right < window.stop:
            right_flank = max(0, min(right_flank, seq_right - tx_end))
        effective_buffer = max(0, min(buffer_bp, left_flank, right_flank))

    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise SequencesError(
            "sequences_timeout",
            "sequences request exceeded deadline",
            status=504,
            **_group_error_extra(group),
        )
    gffread_timeout = min(float(settings.GFFREAD_TIMEOUT_S), max(1.0, remaining))

    try:
        run_gffread(
            gff_path=gff_path,
            fasta_path=fa_path,
            out_path=out_raw,
            action=action,
            buffer_bp=effective_buffer,
            include_exon_coords=include_exon_coords,
            timeout_s=gffread_timeout,
        )
    except GffreadError as exc:
        status = 503 if exc.code == "gffread_unavailable" else 500
        if exc.code == "buffer_not_supported":
            status = 422
        raise SequencesError(
            exc.code,
            str(exc),
            status=status,
            **_merge_extra(_group_error_extra(group), exc.extra),
        ) from exc

    expected = len(group.items)
    got = count_fasta_records(out_raw)
    if got != expected:
        raise SequencesError(
            "gffread_failed",
            f"gffread returned {got} records, expected {expected}",
            status=500,
            **_merge_extra(
                _group_error_extra(group),
                {"expected": expected, "got": got},
            ),
        )

    try:
        by_tid = split_fasta_to_record_files(
            out_raw,
            records_dir,
            header_for_id=header_for_id,
        )
    except GffreadError as exc:
        raise SequencesError(
            exc.code,
            str(exc),
            status=500,
            **_merge_extra(_group_error_extra(group), exc.extra),
        ) from exc

    for item in group.items:
        if item.transcript_id not in by_tid:
            raise SequencesError(
                "gffread_failed",
                f"missing FASTA record for {item.transcript_id}",
                status=500,
                **_item_extra(item),
            )
    return by_tid


def _validate_truncation(
    group: GeneGroup,
    *,
    window: GenomicWindow,
    sequence: str,
    buffer_bp: int,
) -> None:
    """Accept short NCBI results only when the shortfall is buffer flank."""
    expected = window.length
    got = len(sequence)
    if got == expected:
        return
    if got > expected:
        raise SequencesError(
            "ncbi_sequence_truncated",
            f"NCBI returned longer sequence than requested ({got} > {expected})",
            status=502,
            annotation_id=group.annotation_id,
            local_id=group.local_id,
            expected=expected,
            got=got,
        )
    # Sequence covers [window.start, window.start + got - 1]
    covered_stop = window.start + got - 1
    for item in group.items:
        model = group.transcripts[item.transcript_id]
        if model.start < window.start or model.end > covered_stop:
            raise SequencesError(
                "ncbi_sequence_truncated",
                (
                    f"NCBI sequence truncated inside transcript "
                    f"{item.transcript_id} (covered through {covered_stop}, "
                    f"transcript ends at {model.end})"
                ),
                status=502,
                **_merge_extra(
                    _item_extra(item),
                    {"expected": expected, "got": got},
                ),
            )
    # Shortfall must lie only in the 3' buffer (efetch clips at chromosome end).
    # 5' is already clipped to 1 when building the window.
    if buffer_bp <= 0:
        raise SequencesError(
            "ncbi_sequence_truncated",
            "NCBI sequence truncated without buffer flank",
            status=502,
            annotation_id=group.annotation_id,
            local_id=group.local_id,
            expected=expected,
            got=got,
        )


def export_sequences(
    request: SequencesRequest,
    *,
    corpus: sqlite3.Connection,
) -> SequencesResult:
    """
    Fail-first bulk export. Caller owns cleaning ``result.temp_dir``.
    """
    deadline = time.monotonic() + float(settings.SEQUENCES_DEADLINE_S)
    groups = _group_items(request.items)

    shard_cache: dict[str, sqlite3.Connection] = {}
    ann_cache: dict[str, dict] = {}
    files_cache: dict[str, AnnotationFiles] = {}
    accession_cache: dict[tuple[str, str], str] = {}

    try:
        # Phase 2: local resolve (no network).
        for group in groups:
            _deadline_check(deadline)
            _resolve_group(
                group,
                corpus=corpus,
                shard_cache=shard_cache,
                ann_cache=ann_cache,
                files_cache=files_cache,
                accession_cache=accession_cache,
                action=request.action,
                buffer_bp=request.buffer_bp,
            )

        _deadline_check(deadline)

        # Phase 3: coalesce + fetch NCBI windows.
        windows = [g.window for g in groups if g.window is not None]
        coalesced = coalesce_windows(
            windows,  # type: ignore[arg-type]
            gap_bp=int(settings.SEQUENCES_COALESCE_GAP_BP),
            max_merged_bp=min(1_000_000, int(settings.SEQUENCES_MAX_WINDOW_BP)),
        )

        client = get_ncbi_client()
        fetched: dict[tuple[str, int, int], str] = {}

        def _fetch_one(cw):
            slice_ = client.fetch_window(
                GenomicWindow(
                    accession=cw.accession,
                    start=cw.start,
                    stop=cw.stop,
                ),
                allow_truncation=True,
            )
            return (cw.accession, cw.start, cw.stop), slice_.sequence

        # Bound parallelism; share the process rate limiter.
        max_workers = min(
            4,
            max(1, len(coalesced)),
            max(1, int(settings.NCBI_WORKERS)),
        )
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {pool.submit(_fetch_one, cw): cw for cw in coalesced}
            try:
                for fut in as_completed(futures):
                    _deadline_check(deadline)
                    try:
                        key, seq = fut.result()
                        fetched[key] = seq
                    except NcbiError as exc:
                        for pending in futures:
                            pending.cancel()
                        raise SequencesError(
                            exc.code,
                            str(exc),
                            status=exc.status,
                            retry_after=exc.retry_after,
                            **_merge_extra(exc.extra),
                        ) from exc
                    except SequencesError:
                        for pending in futures:
                            pending.cancel()
                        raise
            except SequencesError:
                raise

        # Map each group window to a sequence cut from its coalesced parent.
        group_sequences: dict[str, tuple[GenomicWindow, str]] = {}
        for group in groups:
            assert group.window is not None
            w = group.window
            parent_seq = None
            parent_start = None
            for (acc, start, stop), seq in fetched.items():
                if acc != w.accession:
                    continue
                if start <= w.start and stop >= w.stop:
                    parent_seq = seq
                    parent_start = start
                    break
                # Truncated parent may still cover the transcript core.
                if start <= w.start and (start + len(seq) - 1) >= w.start:
                    parent_seq = seq
                    parent_start = start
                    break
            if parent_seq is None or parent_start is None:
                raise SequencesError(
                    "ncbi_error",
                    "internal error: missing fetched window for gene group",
                    status=502,
                    annotation_id=group.annotation_id,
                    local_id=group.local_id,
                )
            # Desired slice may be truncated at the chromosome end.
            available_stop = parent_start + len(parent_seq) - 1
            slice_stop = min(w.stop, available_stop)
            if slice_stop < w.start:
                raise SequencesError(
                    "ncbi_sequence_truncated",
                    "NCBI sequence does not cover the transcript window",
                    status=502,
                    annotation_id=group.annotation_id,
                    local_id=group.local_id,
                )
            try:
                cut = slice_sequence(
                    parent_seq,
                    fetched_start=parent_start,
                    slice_start=w.start,
                    slice_stop=slice_stop,
                )
            except ValueError as exc:
                raise SequencesError(
                    "ncbi_sequence_truncated",
                    str(exc),
                    status=502,
                    annotation_id=group.annotation_id,
                    local_id=group.local_id,
                ) from exc
            effective = GenomicWindow(
                accession=w.accession,
                start=w.start,
                stop=w.start + len(cut) - 1,
            )
            _validate_truncation(
                group,
                window=w,
                sequence=cut,
                buffer_bp=request.buffer_bp,
            )
            group_sequences[group.group_key] = (effective, cut)

        _deadline_check(deadline)

        # Phase 4: gffread per group (monotonic gids), stream-assemble final FASTA.
        temp = TemporaryDirectory(prefix="annotrieve-sequences-")
        work = Path(temp.name)
        # item.index -> path to a single-record FASTA file
        record_paths: dict[int, Path] = {}
        try:
            for group_index, group in enumerate(groups):
                _deadline_check(deadline)
                window, sequence = group_sequences[group.group_key]
                by_tid = _run_group_gffread(
                    group,
                    sequence=sequence,
                    window=window,
                    action=request.action,
                    buffer_bp=request.buffer_bp,
                    include_exon_coords=request.include_exon_coords,
                    work_dir=work,
                    gid=f"g{group_index}",
                    deadline=deadline,
                )
                for item in group.items:
                    record_paths[item.index] = by_tid[item.transcript_id]

            final_path = work / f"annotrieve-{request.action.value}.fa"
            with open(final_path, "wb") as out_fh:
                for item in request.items:
                    src = record_paths[item.index]
                    with open(src, "rb") as in_fh:
                        shutil.copyfileobj(in_fh, out_fh, length=64 * 1024)

            filename = (
                f"annotrieve-{request.action.value}-{len(request.items)}.fa"
            )
            return SequencesResult(
                fasta_path=final_path,
                sequence_count=len(request.items),
                temp_dir=temp,
                filename=filename,
            )
        except Exception:
            temp.cleanup()
            raise
    finally:
        for conn in shard_cache.values():
            try:
                conn.close()
            except sqlite3.Error:
                pass
