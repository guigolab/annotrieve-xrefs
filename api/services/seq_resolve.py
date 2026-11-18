"""Resolve one gene for /sequences: GFF slice -> transcript candidates -> chosen transcripts.

Local only (no network). The same result feeds ``/sequences/plan`` and the export. A per-gene
problem never raises: it is recorded as ``reason`` so the caller can skip the gene (manifest) or
fail the request, as ``on_error`` says. Only conditions that affect every gene (for example the
annotation directory not being configured) raise :class:`SeqError`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from helpers.annotation_files import (
    AnnotationFiles,
    AnnotationFilesError,
    resolve_annotation_files,
    resolve_seqid_accession,
)
from helpers.gff_slice import TranscriptModel, load_features_from_lines, resolve_transcripts
from helpers.tabix_fetch import FetchError, fetch_region_lines
from helpers.transcript_select import Candidate, Rule, gene_extent, list_candidates, select_transcripts
from helpers.windows import GenomicWindow
from services.seq_genes import Gene
from settings import settings


class SeqError(Exception):
    """A condition that makes the whole request fail (HTTP status + API error code)."""

    def __init__(self, code: str, message: str, *, status: int = 500, retry_after=None,
                 **extra) -> None:
        super().__init__(message)
        self.code = code
        self.status = status
        self.retry_after = retry_after
        self.extra = extra


@dataclass
class Resolved:
    gene: Gene
    reason: str = ""                       # "" = ok, else a skip reason code
    message: str = ""
    rule: str = ""
    n_candidates: int = 0
    candidates: list[Candidate] = field(default_factory=list)
    chosen: list[Candidate] = field(default_factory=list)
    models: dict[str, TranscriptModel] = field(default_factory=dict)
    window: GenomicWindow | None = None

    @property
    def ok(self) -> bool:
        return not self.reason

    def fail(self, reason: str, message: str) -> "Resolved":
        self.reason, self.message = reason, message
        return self


_SKIP_MESSAGES = {
    "no_transcript": "the gene has no transcript with exon/CDS structure",
    "no_cds": "no transcript of the gene has a CDS",
    "filtered_out": "no transcript passed the 'require' filters",
}


class GeneResolver:
    """Caches annotation files and accession lookups across the genes of one request."""

    def __init__(self, action: str, rule: Rule, buffer_bp: int) -> None:
        self.action = action
        self.rule = rule
        self.buffer_bp = buffer_bp
        self._files: dict[int, AnnotationFiles | AnnotationFilesError] = {}
        self._accessions: dict[tuple[int, str], str | AnnotationFilesError] = {}

    def _annotation_files(self, gene: Gene) -> AnnotationFiles:
        cached = self._files.get(gene.ord)
        if cached is None:
            try:
                cached = resolve_annotation_files(
                    taxid=gene.taxid, assembly_accession=gene.accession, gff_path=gene.gff_path)
            except AnnotationFilesError as exc:
                cached = exc
            self._files[gene.ord] = cached
        if isinstance(cached, AnnotationFilesError):
            raise cached
        return cached

    def _accession(self, files: AnnotationFiles, gene: Gene) -> str:
        key = (gene.ord, gene.seqid)
        cached = self._accessions.get(key)
        if cached is None:
            try:
                cached = resolve_seqid_accession(files, gene.seqid)
            except AnnotationFilesError as exc:
                cached = exc
            self._accessions[key] = cached
        if isinstance(cached, AnnotationFilesError):
            raise cached
        return cached

    def resolve(self, gene: Gene) -> Resolved:
        res = Resolved(gene, rule="pinned" if gene.pinned else self.rule.name)
        if not gene.gene_id:
            return res.fail("gene_id_missing", "the locus has no source gene id to look up in the GFF")
        try:
            files = self._annotation_files(gene)
        except AnnotationFilesError as exc:
            if exc.code == "annotations_unavailable":
                raise SeqError("annotations_unavailable", str(exc), status=503) from exc
            return res.fail("gff_unavailable", str(exc))

        try:
            region = (gene.start, gene.end)
            feats = load_features_from_lines(
                fetch_region_lines(files, seqid=gene.seqid, start=region[0], end=region[1]))
            ext = gene_extent(feats, gene.gene_id)
            if ext and (ext[0] < region[0] or ext[1] > region[1]):
                # A transcript reaches beyond the gene: fetch again so its exons are complete.
                region = (min(region[0], ext[0]), max(region[1], ext[1]))
                feats = load_features_from_lines(
                    fetch_region_lines(files, seqid=gene.seqid, start=region[0], end=region[1]))
        except FetchError as exc:
            reason = "seqid_not_resolved" if exc.code == "seqid_not_resolved" else "gff_unavailable"
            return res.fail(reason, str(exc))

        cands = list_candidates(feats, gene.gene_id)
        res.candidates, res.n_candidates = cands, len(cands)
        if gene.pinned:
            match = [c for c in cands if c.transcript_id == gene.pinned]
            if not match:
                return res.fail("transcript_not_found",
                                f"transcript {gene.pinned!r} not found in the gene's GFF slice")
            if self.action in ("cds", "protein") and not match[0].coding:
                return res.fail("no_cds", f"transcript {gene.pinned!r} has no CDS")
            res.chosen = match
        else:
            sel = select_transcripts(cands, self.rule, self.action)
            if sel.skipped_reason:
                return res.fail(sel.skipped_reason, _SKIP_MESSAGES.get(sel.skipped_reason, ""))
            res.chosen = sel.chosen

        try:
            res.models = resolve_transcripts(
                feats, source_gene_id=gene.gene_id, transcript_ids=[c.transcript_id for c in res.chosen])
        except KeyError as exc:
            return res.fail("transcript_not_found", f"transcript {exc.args[0]!r} not found")
        for tid, model in res.models.items():
            if (model.strand or "").strip() not in ("+", "-"):
                return res.fail("transcript_strand_unknown", f"transcript {tid!r} has no usable strand")
            if not model.has_structure:
                return res.fail("no_transcript", f"transcript {tid!r} has no exon/CDS lines")

        try:
            accession = self._accession(files, gene)
        except AnnotationFilesError as exc:
            return res.fail("seqid_not_resolved", str(exc))
        start = max(1, min(m.start for m in res.models.values()) - self.buffer_bp)
        stop = max(m.end for m in res.models.values()) + self.buffer_bp
        if stop - start + 1 > int(settings.SEQUENCES_MAX_WINDOW_BP):
            return res.fail("window_too_large",
                            f"fetch window {stop - start + 1} bp exceeds {settings.SEQUENCES_MAX_WINDOW_BP}")
        res.window = GenomicWindow(accession=accession, start=start, stop=stop)
        return res
