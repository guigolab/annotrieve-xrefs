"""Batch export for /sequences: resolve genes, fetch genomic windows from NCBI, run gffread.

Genes are processed in small batches (``SEQUENCES_BATCH_GENES``), so memory is bounded by one
batch of NCBI windows and a request of thousands of genes reports progress. Output order is the
order of the input genes; within a gene, the order of the selected transcripts.

A gene that cannot be processed is *skipped* (``on_error=skip``): it gets one manifest row with a
reason code and the export continues. With ``on_error=fail`` the first such gene aborts the
request with a :class:`SeqError`. Conditions that hit every gene (no gffread binary, deadline) are
always fatal.
"""
from __future__ import annotations

import logging
import os
import shutil
import tempfile
import time
from collections import Counter
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from helpers.gff_slice import dedupe_gff_lines, shift_gff_lines, write_fasta, write_gff
from helpers.gffread_runner import (
    GffreadError,
    SequenceAction,
    count_fasta_records,
    run_gffread,
    split_fasta_to_record_files,
)
from helpers.ncbi_eutils import NcbiError, get_ncbi_client
from helpers.windows import CoalescedWindow, GenomicWindow, coalesce_windows, slice_sequence
from services.seq_genes import Gene
from services.seq_request import SeqRequest
from services.seq_resolve import GeneResolver, Resolved, SeqError
from settings import settings

logger = logging.getLogger(__name__)

MANIFEST_COLUMNS = ("locus_ref", "annotation_id", "symbol", "status", "transcript_id", "rule",
                    "n_candidates", "reason", "seq_length", "message")
_STATUS_BY_REASON = {
    "transcript_not_found": 404, "gff_unavailable": 404, "annotation_not_found": 404,
    "locus_not_found": 404, "ncbi_error": 502, "ncbi_sequence_truncated": 502,
}

# (locus_ref, gene or None, reason if the gene could not even be looked up)
Item = tuple[str, "Gene | None", str]
Progress = Callable[[int, "int | None", Counter], None]


@dataclass
class ExportResult:
    fasta_path: Path
    manifest_path: Path
    n_sequences: int
    n_genes: int
    n_skipped: int
    reasons: dict[str, int]


def _clean(text: str) -> str:
    out = "".join("_" if (c.isspace() or c == "|") else c for c in (text or ""))
    return out or "NA"


def _defline(res: Resolved, tid: str, buffer_bp: int) -> str:
    g, model = res.gene, res.models[tid]
    tail = f" {g.seqid}:{model.start}-{model.end}({model.strand})"
    if buffer_bp:
        tail += f" buffer={buffer_bp}"
    return "|".join((_clean(g.symbol or g.gene_id), _clean(g.organism), _clean(g.accession),
                     _clean(g.annotation_id), _clean(tid))) + tail


def _fatal(reason: str, message: str, locus_ref: str, **extra) -> SeqError:
    return SeqError(reason, message or reason, status=_STATUS_BY_REASON.get(reason, 422),
                    locus_ref=locus_ref, **extra)


def _fetch_windows(windows: list[GenomicWindow], deadline: float
                   ) -> tuple[dict[GenomicWindow, CoalescedWindow],
                              dict[tuple[str, int, int], "str | NcbiError"]]:
    coalesced = coalesce_windows(
        windows, gap_bp=int(settings.SEQUENCES_COALESCE_GAP_BP),
        max_merged_bp=min(1_000_000, int(settings.SEQUENCES_MAX_WINDOW_BP)))
    parent = {m: cw for cw in coalesced for m in cw.members}
    client = get_ncbi_client()

    def one(cw: CoalescedWindow):
        key = (cw.accession, cw.start, cw.stop)
        try:
            return key, client.fetch_window(
                GenomicWindow(cw.accession, cw.start, cw.stop), allow_truncation=True).sequence
        except NcbiError as exc:
            return key, exc

    fetched: dict[tuple[str, int, int], str | NcbiError] = {}
    workers = min(4, max(1, len(coalesced)), max(1, int(settings.NCBI_WORKERS)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for key, value in pool.map(one, coalesced):
            if time.monotonic() > deadline:
                raise SeqError("sequences_timeout", "sequences request exceeded its deadline",
                               status=504)
            fetched[key] = value
    return parent, fetched


def _cut(res: Resolved, parent: CoalescedWindow, fetched, buffer_bp: int
         ) -> tuple[GenomicWindow, str] | tuple[str, str]:
    """(window, sequence) for one gene, or (reason, message) when NCBI data is unusable."""
    w = res.window
    value = fetched[(parent.accession, parent.start, parent.stop)]
    if isinstance(value, NcbiError):
        return "ncbi_error", str(value)
    available_stop = parent.start + len(value) - 1
    slice_stop = min(w.stop, available_stop)
    if slice_stop < w.start:
        return "ncbi_sequence_truncated", "NCBI sequence does not cover the transcript window"
    cut = slice_sequence(value, fetched_start=parent.start, slice_start=w.start,
                         slice_stop=slice_stop)
    if len(cut) < w.length:
        covered_stop = w.start + len(cut) - 1
        for tid, model in res.models.items():
            if model.start < w.start or model.end > covered_stop:
                return "ncbi_sequence_truncated", (
                    f"NCBI sequence truncated inside transcript {tid} (ends at {covered_stop})")
        if buffer_bp <= 0:
            return "ncbi_sequence_truncated", "NCBI sequence truncated without buffer flank"
    return GenomicWindow(w.accession, w.start, w.start + len(cut) - 1), cut


def _run_gffread(res: Resolved, window: GenomicWindow, sequence: str, req: SeqRequest,
                 work: Path, gid: str, deadline: float) -> dict[str, Path]:
    chosen_ids = [c.transcript_id for c in res.chosen]
    lines: list[str] = []
    for tid in chosen_ids:
        lines.extend(res.models[tid].lines)
    shifted = shift_gff_lines(dedupe_gff_lines(lines), new_seqid=gid, origin=window.start)
    gff, fa, raw = work / f"{gid}.gff3", work / f"{gid}.fa", work / f"{gid}.raw.fa"
    write_gff(gff, shifted)
    write_fasta(fa, seqid=gid, sequence=sequence)

    buffer_bp = req.buffer_bp
    if buffer_bp and req.action is SequenceAction.transcript:
        tx_start = min(m.start for m in res.models.values())
        tx_end = max(m.end for m in res.models.values())
        right = window.stop - tx_end
        buffer_bp = max(0, min(buffer_bp, tx_start - window.start, right))
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise SeqError("sequences_timeout", "sequences request exceeded its deadline", status=504)
    run_gffread(gff_path=gff, fasta_path=fa, out_path=raw, action=req.action,
                buffer_bp=buffer_bp, include_exon_coords=req.include_exon_coords,
                timeout_s=min(float(settings.GFFREAD_TIMEOUT_S), max(1.0, remaining)))
    got = count_fasta_records(raw)
    if got != len(chosen_ids):
        raise GffreadError("gffread_failed", f"gffread returned {got} records, expected "
                           f"{len(chosen_ids)}")
    headers = {tid: _defline(res, tid, req.buffer_bp) for tid in chosen_ids}
    return split_fasta_to_record_files(raw, work / f"{gid}.records", header_for_id=headers)


def _row(*cols) -> str:
    return "\t".join(str(c).replace("\t", " ").replace("\n", " ") for c in cols) + "\n"


def run_export(req: SeqRequest, items: Iterable[Item], *, out_dir: Path,
               total: int | None = None, progress: Progress | None = None,
               deadline: float | None = None) -> ExportResult:
    deadline = deadline if deadline is not None else (
        time.monotonic() + float(settings.SEQUENCES_DEADLINE_S))
    out_dir.mkdir(parents=True, exist_ok=True)
    resolver = GeneResolver(req.action.value, req.rule, req.buffer_bp)
    fasta_part, manifest_part = out_dir / "result.fa.part", out_dir / "manifest.tsv.part"
    counts: Counter = Counter()
    reasons: Counter = Counter()
    n_genes = n_seq = 0
    gene_counter = 0
    work_root = Path(tempfile.mkdtemp(prefix="seq-", dir=out_dir))

    def skip(manifest, ref: str, gene: Gene | None, reason: str, message: str,
             res: Resolved | None = None) -> None:
        if req.on_error == "fail":
            raise _fatal(reason, message, ref)
        counts["skipped"] += 1
        reasons[reason] += 1
        manifest.write(_row(ref, gene.annotation_id if gene else ref.rpartition(":")[0],
                            gene.symbol if gene else "", "skipped", "",
                            res.rule if res else "", res.n_candidates if res else 0,
                            reason, 0, message))

    try:
        with open(fasta_part, "wb") as fasta, open(manifest_part, "w", encoding="utf-8") as manifest:
            manifest.write("\t".join(MANIFEST_COLUMNS) + "\n")

            def flush(batch: list[Item]) -> None:
                nonlocal n_seq, gene_counter
                if time.monotonic() > deadline:
                    raise SeqError("sequences_timeout", "sequences request exceeded its deadline",
                                   status=504)
                resolved: list[tuple[str, Gene | None, Resolved | None, str, str]] = []
                for ref, gene, reason in batch:
                    if gene is None:
                        resolved.append((ref, None, None, reason, reason.replace("_", " ")))
                        continue
                    res = resolver.resolve(gene)
                    resolved.append((ref, gene, res, res.reason, res.message))

                ok = [r for _, _, r, reason, _ in resolved if r is not None and not reason]
                parent, fetched = ({}, {})
                if ok:
                    parent, fetched = _fetch_windows([r.window for r in ok], deadline)

                for ref, gene, res, reason, message in resolved:
                    if reason or res is None:
                        skip(manifest, ref, gene, reason, message, res)
                        continue
                    cut = _cut(res, parent[res.window], fetched, req.buffer_bp)
                    if isinstance(cut[0], str):
                        skip(manifest, ref, gene, cut[0], cut[1], res)
                        continue
                    window, sequence = cut
                    gene_counter += 1
                    try:
                        records = _run_gffread(res, window, sequence, req, work_root,
                                               f"g{gene_counter}", deadline)
                    except GffreadError as exc:
                        if exc.code == "gffread_unavailable":
                            raise SeqError(exc.code, str(exc), status=503) from exc
                        skip(manifest, ref, gene, "gffread_failed", str(exc), res)
                        continue
                    counts["ok"] += 1
                    for cand in res.chosen:
                        path = records[cand.transcript_id]
                        length = 0
                        with open(path, "rb") as fh:
                            for line in fh:
                                fasta.write(line)
                                if not line.startswith(b">"):
                                    length += len(line.strip())
                        n_seq += 1
                        manifest.write(_row(ref, gene.annotation_id, gene.symbol, "ok",
                                            cand.transcript_id, res.rule, res.n_candidates, "",
                                            length, ""))
                    shutil.rmtree(work_root / f"g{gene_counter}.records", ignore_errors=True)

            batch: list[Item] = []
            for item in items:
                batch.append(item)
                if len(batch) >= int(settings.SEQUENCES_BATCH_GENES):
                    n_genes += len(batch)
                    flush(batch)
                    batch = []
                    if progress:
                        progress(n_genes, total, counts)
            if batch:
                n_genes += len(batch)
                flush(batch)
                if progress:
                    progress(n_genes, total, counts)
    except BaseException:
        fasta_part.unlink(missing_ok=True)
        manifest_part.unlink(missing_ok=True)
        raise
    finally:
        shutil.rmtree(work_root, ignore_errors=True)

    fasta_path, manifest_path = out_dir / "result.fa", out_dir / "manifest.tsv"
    os.replace(fasta_part, fasta_path)
    os.replace(manifest_part, manifest_path)
    return ExportResult(fasta_path, manifest_path, n_seq, n_genes, counts["skipped"],
                        dict(reasons))
