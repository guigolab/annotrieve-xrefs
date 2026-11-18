"""Choose which transcripts of a gene a /sequences request extracts.

``transcript=`` is a preset name or a rule object::

    "longest"
    {"rank_by": ["coding:desc", "cds_len:desc"], "require": ["coding", "tag:basic"], "keep": 1}

Rank keys (each ``key:asc`` or ``key:desc``; a bare key means desc):

    coding        has CDS (1/0)             spliced_len  merged exon/CDS/UTR length
    cds_len       summed CDS length         span         genomic start..end
    n_exons       exon count (CDS segments if the transcript has no exon lines)
    utr_len       spliced_len - cds_len

Filters (``require``): ``coding``, ``tag:<value>`` (GFF ``tag`` attribute), ``biotype:<value>``.

Ties always go to the lowest transcript id, so a request is reproducible. A transcript is
*eligible* when it can give the action's sequence: CDS for ``cds``/``protein``, any structure for
``transcript``. Eligibility is applied before ranking, so a non-coding transcript never wins a
``cds`` request.

A gene whose structural features hang directly on the gene (no mRNA) is represented by one
synthetic transcript whose id is the gene id, as in ``gff_slice.resolve_transcripts``.

There is deliberately no "canonical" key: the flag exists in one data source only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from helpers.gff_slice import GffFeature, _is_gene, _is_structural

RANK_KEYS = ("coding", "cds_len", "spliced_len", "span", "n_exons", "utr_len")
PRESETS = ("longest", "longest_spliced", "longest_cds", "shortest_cds", "most_exons",
           "fewest_exons", "all")
MAX_KEEP = 50


class RuleError(ValueError):
    """Invalid ``transcript`` parameter. ``code`` is an API error code."""

    def __init__(self, code: str, message: str, **extra) -> None:
        super().__init__(message)
        self.code = code
        self.extra = extra


@dataclass(frozen=True)
class Candidate:
    transcript_id: str
    synthetic: bool
    spliced_len: int
    cds_len: int
    n_exons: int
    span: int
    biotype: str
    tags: tuple[str, ...]

    @property
    def coding(self) -> bool:
        return self.cds_len > 0

    @property
    def utr_len(self) -> int:
        return max(self.spliced_len - self.cds_len, 0)

    def metrics(self) -> dict:
        return {"coding": self.coding, "cds_len": self.cds_len, "spliced_len": self.spliced_len,
                "span": self.span, "n_exons": self.n_exons, "utr_len": self.utr_len,
                "biotype": self.biotype, "tags": list(self.tags)}


@dataclass(frozen=True)
class Rule:
    name: str                                   # preset name or "custom"
    rank_by: tuple[tuple[str, str], ...]        # (key, "asc"|"desc")
    require: tuple[str, ...] = ()
    keep: int | None = 1                        # None = all


@dataclass(frozen=True)
class Selection:
    chosen: list[Candidate]
    skipped_reason: str = ""                    # no_transcript | no_cds | filtered_out
    rule: str = ""
    n_candidates: int = 0
    ranked: list[Candidate] = field(default_factory=list)   # eligible, best first (for /sequence-plans)


def normalise_tag(value: str) -> str:
    return value.strip().casefold().replace(" ", "_").replace("-", "_")


def feature_tags(feat: GffFeature) -> tuple[str, ...]:
    """Normalised ``tag`` values of one feature (attribute key matched case-insensitively)."""
    out: list[str] = []
    for key, value in feat.attrs.items():
        if key.casefold() in ("tag", "tags"):
            out.extend(normalise_tag(part) for part in value.split(",") if part.strip())
    return tuple(out)


def _biotype(feat: GffFeature) -> str:
    for key in ("biotype", "transcript_biotype", "gene_biotype"):
        for k, v in feat.attrs.items():
            if k.casefold() == key and v:
                return v.strip().casefold()
    return feat.type.casefold()


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


def _measure(parts: list[GffFeature]) -> tuple[int, int, int]:
    """(spliced_len, cds_len, n_exons) of a transcript's structural children."""
    body = [c for c in parts if c.type.casefold() != "intron"]
    cds = [c for c in body if c.type.casefold() == "cds"]
    exons = [c for c in body if c.type.casefold() == "exon"]
    return (_merged_length([(c.start, c.end) for c in body]),
            _merged_length([(c.start, c.end) for c in cds]),
            len(exons) or len(cds))


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
        spliced, cds, n_exons = _measure(parts)
        out.append(Candidate(feat.feature_id, False, spliced, cds, n_exons,
                             feat.end - feat.start + 1, _biotype(feat), feature_tags(feat)))
    if not out:
        direct = [c for c in children.get(source_gene_id, []) if _is_structural(c.type)]
        if direct:
            spliced, cds, n_exons = _measure(direct)
            span = max(c.end for c in direct) - min(c.start for c in direct) + 1
            out.append(Candidate(source_gene_id, True, spliced, cds, n_exons, span, "", ()))
    return sorted(out, key=lambda c: c.transcript_id)


def gene_extent(features: list[GffFeature], source_gene_id: str) -> tuple[int, int] | None:
    """(start, end) covering the gene's own transcript lines, or None if it has none.

    A tabix slice returns only features that overlap the queried region, so a transcript whose
    exons lie partly outside the gene bounds would come back incomplete. The caller compares this
    extent with the region it fetched and fetches again when it is larger.
    """
    spans = [(f.start, f.end) for f in features
             if source_gene_id in f.parents and f.feature_id
             and not _is_structural(f.type) and not _is_gene(f.type)]
    if not spans:
        return None
    return min(s for s, _ in spans), max(e for _, e in spans)


# ---------------------------------------------------------------- rule parsing

def _preset(name: str, action: str) -> Rule:
    coding = action in ("cds", "protein")
    table = {
        "longest": ((("cds_len", "desc"),) if coding
                    else (("coding", "desc"), ("spliced_len", "desc"))),
        "longest_spliced": (("spliced_len", "desc"),),
        "longest_cds": (("cds_len", "desc"),),
        "shortest_cds": (("cds_len", "asc"),),
        "most_exons": (("n_exons", "desc"),),
        "fewest_exons": (("n_exons", "asc"),),
        "all": (),
    }
    require = ("coding",) if name in ("longest_cds", "shortest_cds") else ()
    return Rule(name, table[name], require, None if name == "all" else 1)


def _parse_rank(item: object) -> tuple[str, str]:
    if not isinstance(item, str):
        raise RuleError("invalid_transcript_rule", "rank_by entries must be strings like 'cds_len:desc'")
    key, _, direction = item.partition(":")
    key, direction = key.strip(), (direction.strip() or "desc")
    if key not in RANK_KEYS:
        raise RuleError("unknown_rank_key", f"unknown rank key {key!r}", allowed=list(RANK_KEYS))
    if direction not in ("asc", "desc"):
        raise RuleError("invalid_transcript_rule", f"direction must be asc or desc, got {direction!r}")
    return key, direction


def _check_require(item: object) -> str:
    if not isinstance(item, str):
        raise RuleError("invalid_transcript_rule", "require entries must be strings")
    text = item.strip()
    if text == "coding":
        return text
    kind, _, value = text.partition(":")
    if kind == "tag" and value.strip():
        return "tag:" + normalise_tag(value)
    if kind == "biotype" and value.strip():
        return "biotype:" + value.strip().casefold()
    raise RuleError("unknown_require", f"unknown require filter {text!r}",
                    allowed=["coding", "tag:<value>", "biotype:<value>"])


def parse_rule(spec: object, action: str) -> Rule:
    """``spec``: None (default), a preset name, or a dict with rank_by/require/keep."""
    if spec is None:
        return _preset("longest", action)
    if isinstance(spec, str):
        if spec not in PRESETS:
            raise RuleError("unknown_transcript_preset", f"unknown transcript preset {spec!r}",
                            allowed=list(PRESETS))
        return _preset(spec, action)
    if not isinstance(spec, dict):
        raise RuleError("invalid_transcript_rule", "transcript must be a preset name or an object")
    unknown = set(spec) - {"rank_by", "require", "keep"}
    if unknown:
        raise RuleError("invalid_transcript_rule", f"unknown transcript rule fields {sorted(unknown)}")
    rank = tuple(_parse_rank(i) for i in spec.get("rank_by") or [])
    require = tuple(_check_require(i) for i in spec.get("require") or [])
    keep = spec.get("keep", 1)
    if keep == "all":
        keep = None
    elif not isinstance(keep, int) or isinstance(keep, bool) or not 1 <= keep <= MAX_KEEP:
        raise RuleError("invalid_transcript_rule", f"keep must be 'all' or an integer 1..{MAX_KEEP}")
    return Rule("custom", rank, require, keep)


# ------------------------------------------------------------------- selection

def _value(c: Candidate, key: str) -> int:
    return int(getattr(c, key))


def _passes(c: Candidate, req: str) -> bool:
    if req == "coding":
        return c.coding
    kind, _, value = req.partition(":")
    if kind == "tag":
        return value in c.tags
    return c.biotype == value


def select_transcripts(cands: list[Candidate], rule: Rule, action: str) -> Selection:
    """Apply ``rule`` to the candidates of one gene for ``action`` (cds | protein | transcript)."""
    n = len(cands)
    if not cands:
        return Selection([], "no_transcript", rule.name, 0)
    need_cds = action in ("cds", "protein")
    eligible = [c for c in cands if (c.coding if need_cds else c.spliced_len > 0)]
    if not eligible:
        return Selection([], "no_cds" if need_cds else "no_transcript", rule.name, n)
    kept = [c for c in eligible if all(_passes(c, r) for r in rule.require)]
    if not kept:
        return Selection([], "filtered_out", rule.name, n)

    def sort_key(c: Candidate):
        parts = []
        for key, direction in rule.rank_by:
            v = _value(c, key)
            parts.append(-v if direction == "desc" else v)
        return (*parts, c.transcript_id)

    ranked = sorted(kept, key=sort_key)
    chosen = ranked if rule.keep is None else ranked[: rule.keep]
    return Selection(chosen, "", rule.name, n, ranked)
