"""Choose which transcripts of a gene a /sequences request extracts.

Modes (``transcript=``):

- ``all``                 every eligible transcript of the gene.
- ``longest``             one per gene: the longest by the metric of the action.
- ``canonical_or_longest`` one per gene: the canonical transcript, else the longest.
- ``canonical_only``      one per gene: the canonical transcript; genes without one are dropped.

Length metric (the same for every mode that compares transcripts):

- ``cds`` / ``protein`` -> summed CDS length (transcripts without CDS are not eligible);
- ``transcript``        -> summed length of the merged exon/CDS/UTR intervals (spliced length).

Ties are broken by the lowest transcript id, so the result is deterministic. When several
transcripts of a gene are flagged canonical (for example MANE Select and Ensembl canonical on
different transcripts), the longest of the flagged ones is chosen.

A gene whose structural features hang directly on the gene (no mRNA) is represented by one
synthetic transcript whose id is the gene id, as in ``gff_slice.resolve_transcripts``.

Canonical flags are read from the ``tag`` attribute of the transcript line. The recognised values
are listed in ``CANONICAL_TAGS`` after normalisation (lower case, spaces and dashes as ``_``).
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from helpers.gff_slice import GffFeature, _is_gene, _is_structural

CANONICAL_TAGS = frozenset({"ensembl_canonical", "mane_select", "refseq_select"})


class TranscriptMode(str, Enum):
    all = "all"
    longest = "longest"
    canonical_or_longest = "canonical_or_longest"
    canonical_only = "canonical_only"


@dataclass(frozen=True)
class Candidate:
    transcript_id: str
    synthetic: bool
    spliced_len: int
    cds_len: int
    canonical: bool
    tags: tuple[str, ...]


@dataclass(frozen=True)
class Selection:
    chosen: list[Candidate]
    # Why nothing was chosen ("no_transcript", "no_cds", "no_canonical") or "" when something was.
    skipped_reason: str = ""
    # How the single transcript was picked: "canonical", "longest", "all" ("" when skipped).
    rule: str = ""


def normalise_tag(value: str) -> str:
    return value.strip().casefold().replace(" ", "_").replace("-", "_")


def feature_tags(feat: GffFeature) -> tuple[str, ...]:
    """Normalised ``tag`` values of one feature (attribute key matched case-insensitively)."""
    out: list[str] = []
    for key, value in feat.attrs.items():
        if key.casefold() not in ("tag", "tags"):
            continue
        out.extend(normalise_tag(part) for part in value.split(",") if part.strip())
    return tuple(out)


def _merged_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    cur_s = cur_e = None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def _lengths(children: list[GffFeature]) -> tuple[int, int]:
    structural = [(c.start, c.end) for c in children if _is_structural(c.type)
                  and c.type.casefold() != "intron"]
    cds = [(c.start, c.end) for c in children if c.type.casefold() == "cds"]
    return _merged_length(structural), _merged_length(cds)


def list_candidates(features: list[GffFeature], source_gene_id: str) -> list[Candidate]:
    """All transcripts of ``source_gene_id`` found in ``features`` (a tabix slice of the gene)."""
    children: dict[str, list[GffFeature]] = {}
    for feat in features:
        for parent in feat.parents:
            children.setdefault(parent, []).append(feat)

    out: list[Candidate] = []
    for feat in children.get(source_gene_id, []):
        if not feat.feature_id or _is_structural(feat.type) or _is_gene(feat.type):
            continue
        parts = [c for c in children.get(feat.feature_id, []) if _is_structural(c.type)]
        if not parts:
            continue
        spliced, cds = _lengths(parts)
        tags = feature_tags(feat)
        out.append(Candidate(feat.feature_id, False, spliced, cds,
                             any(t in CANONICAL_TAGS for t in tags), tags))
    if not out:
        direct = [c for c in children.get(source_gene_id, []) if _is_structural(c.type)]
        if direct:
            spliced, cds = _lengths(direct)
            out.append(Candidate(source_gene_id, True, spliced, cds, False, ()))
    return sorted(out, key=lambda c: c.transcript_id)


def _metric(c: Candidate, action: str) -> int:
    return c.cds_len if action in ("cds", "protein") else c.spliced_len


def _longest(cands: list[Candidate], action: str) -> Candidate:
    return min(cands, key=lambda c: (-_metric(c, action), c.transcript_id))


def select_transcripts(cands: list[Candidate], mode: TranscriptMode | str, action: str) -> Selection:
    """Apply ``mode`` to the candidates of one gene for ``action`` (cds | protein | transcript)."""
    mode = TranscriptMode(mode)
    if not cands:
        return Selection([], "no_transcript")
    eligible = [c for c in cands if _metric(c, action) > 0]
    if not eligible:
        return Selection([], "no_cds" if action in ("cds", "protein") else "no_transcript")
    if mode is TranscriptMode.all:
        return Selection(eligible, "", "all")
    canonical = [c for c in eligible if c.canonical]
    if mode is TranscriptMode.longest:
        return Selection([_longest(eligible, action)], "", "longest")
    if canonical:
        return Selection([_longest(canonical, action)], "", "canonical")
    if mode is TranscriptMode.canonical_or_longest:
        return Selection([_longest(eligible, action)], "", "longest")
    return Selection([], "no_canonical")
