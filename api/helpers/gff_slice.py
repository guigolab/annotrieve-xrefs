"""Parse a tabix GFF gene-window and extract transcript subtrees for gffread."""
from __future__ import annotations

from dataclasses import dataclass, field
from urllib.parse import unquote

# Structural children of a transcript / gene (front/lib/genes/gff.ts).
_STRUCTURAL = frozenset(
    {
        "exon",
        "cds",
        "utr",
        "five_prime_utr",
        "three_prime_utr",
        "5utr",
        "3utr",
        "intron",
    }
)
_GENE_TYPES = frozenset({"gene", "pseudogene"})


def _is_structural(ftype: str) -> bool:
    t = ftype.casefold()
    if t in _STRUCTURAL:
        return True
    return t.endswith("_utr")


def _is_gene(ftype: str) -> bool:
    return ftype.casefold() in _GENE_TYPES


@dataclass
class GffFeature:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes_raw: str
    attrs: dict[str, str]
    feature_id: str
    parents: list[str]
    raw_line: str


@dataclass
class TranscriptModel:
    transcript_id: str
    seqid: str
    start: int
    end: int
    strand: str
    synthetic: bool
    has_cds: bool
    has_structure: bool
    # Lines to write for gffread (transcript + children), already Parent-fixed.
    lines: list[str] = field(default_factory=list)


def parse_gff_attributes(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in (raw or "").split(";"):
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        key = key.strip()
        if not key:
            continue
        value = value.strip()
        try:
            value = unquote(value)
        except Exception:
            pass
        out[key] = value
    return out


def parse_parents(attrs: dict[str, str]) -> list[str]:
    raw = attrs.get("Parent") or attrs.get("parent") or ""
    return [p.strip() for p in raw.split(",") if p.strip()]


def parse_gff_line(line: str) -> GffFeature | None:
    if not line or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9:
        return None
    attrs = parse_gff_attributes(parts[8])
    fid = attrs.get("ID") or attrs.get("Id") or attrs.get("id") or ""
    try:
        start = int(parts[3])
        end = int(parts[4])
    except ValueError:
        return None
    return GffFeature(
        seqid=parts[0],
        source=parts[1],
        type=parts[2],
        start=start,
        end=end,
        score=parts[5],
        strand=parts[6],
        phase=parts[7],
        attributes_raw=parts[8],
        attrs=attrs,
        feature_id=fid,
        parents=parse_parents(attrs),
        raw_line=line.rstrip("\n"),
    )


def _rewrite_attributes(
    attrs: dict[str, str],
    *,
    set_id: str | None = None,
    set_parents: list[str] | None = None,
) -> str:
    ordered: list[tuple[str, str]] = []
    seen: set[str] = set()
    if set_id is not None:
        ordered.append(("ID", set_id))
        seen.add("ID")
    if set_parents is not None:
        ordered.append(("Parent", ",".join(set_parents)))
        seen.add("Parent")
    for key, value in attrs.items():
        if key in seen or key.casefold() in {s.casefold() for s in seen}:
            continue
        if set_id is not None and key.casefold() == "id":
            continue
        if set_parents is not None and key.casefold() == "parent":
            continue
        ordered.append((key, value))
        seen.add(key)
    return ";".join(f"{k}={v}" for k, v in ordered)


def _format_line(
    feat: GffFeature,
    *,
    seqid: str,
    start: int,
    end: int,
    ftype: str | None = None,
    attrs_raw: str | None = None,
    strand: str | None = None,
) -> str:
    return "\t".join(
        [
            seqid,
            feat.source,
            ftype if ftype is not None else feat.type,
            str(start),
            str(end),
            feat.score,
            strand if strand is not None else feat.strand,
            feat.phase,
            attrs_raw if attrs_raw is not None else feat.attributes_raw,
        ]
    )


def load_features_from_lines(lines: list[str]) -> list[GffFeature]:
    features: list[GffFeature] = []
    for line in lines:
        feat = parse_gff_line(line)
        if feat is not None:
            features.append(feat)
    return features


def resolve_transcripts(
    features: list[GffFeature],
    *,
    source_gene_id: str,
    transcript_ids: list[str],
) -> dict[str, TranscriptModel]:
    """
    Build transcript models for each requested id from an in-memory feature list.

    Synthetic: when ``transcript_id == source_gene_id``, structural children of
    the gene are re-parented under a synthesized mRNA with that id.
    """
    by_id: dict[str, GffFeature] = {}
    children: dict[str, list[GffFeature]] = {}
    for feat in features:
        if feat.feature_id:
            by_id[feat.feature_id] = feat
        for parent in feat.parents:
            children.setdefault(parent, []).append(feat)

    out: dict[str, TranscriptModel] = {}
    for tid in transcript_ids:
        if tid in out:
            continue
        if tid == source_gene_id:
            model = _synthetic_transcript(
                features,
                by_id=by_id,
                children=children,
                source_gene_id=source_gene_id,
            )
        else:
            model = _normal_transcript(
                by_id=by_id,
                children=children,
                transcript_id=tid,
            )
        out[tid] = model
    return out


def _has_cds(feats: list[GffFeature]) -> bool:
    return any(f.type.casefold() == "cds" for f in feats)


def _normal_transcript(
    *,
    by_id: dict[str, GffFeature],
    children: dict[str, list[GffFeature]],
    transcript_id: str,
) -> TranscriptModel:
    tx = by_id.get(transcript_id)
    if tx is None or _is_structural(tx.type) or _is_gene(tx.type):
        raise KeyError(transcript_id)
    parts = [
        c
        for c in children.get(transcript_id, [])
        if transcript_id in c.parents and _is_structural(c.type)
    ]
    # Also keep non-structural children that claim this parent (rare).
    extras = [
        c
        for c in children.get(transcript_id, [])
        if transcript_id in c.parents and not _is_structural(c.type) and c is not tx
    ]
    strand = (tx.strand or "").strip()
    lines = [tx.raw_line, *[p.raw_line for p in parts], *[e.raw_line for e in extras]]
    return TranscriptModel(
        transcript_id=transcript_id,
        seqid=tx.seqid,
        start=tx.start,
        end=tx.end,
        strand=strand,
        synthetic=False,
        has_cds=_has_cds(parts),
        has_structure=bool(parts),
        lines=lines,
    )


def _synthetic_transcript(
    features: list[GffFeature],
    *,
    by_id: dict[str, GffFeature],
    children: dict[str, list[GffFeature]],
    source_gene_id: str,
) -> TranscriptModel:
    gene = by_id.get(source_gene_id)
    direct = [
        c
        for c in children.get(source_gene_id, [])
        if source_gene_id in c.parents and _is_structural(c.type)
    ]
    if not direct:
        raise KeyError(source_gene_id)
    start = min(f.start for f in direct)
    end = max(f.end for f in direct)
    strand = (gene.strand if gene else direct[0].strand) or direct[0].strand
    seqid = (gene.seqid if gene else direct[0].seqid)
    # Synthesize an mRNA line so gffread sees a normal transcript.
    # Gene feature itself is omitted from the output; Parent points at it
    # when present so the hierarchy stays GFF3-valid.
    if gene is not None:
        attrs = _rewrite_attributes(
            {"ID": source_gene_id},
            set_id=source_gene_id,
            set_parents=[source_gene_id],
        )
    else:
        attrs = f"ID={source_gene_id}"
    template = gene or direct[0]
    mrna_line = _format_line(
        template,
        seqid=seqid,
        start=start,
        end=end,
        ftype="mRNA",
        attrs_raw=attrs,
        strand=strand,
    )
    child_lines: list[str] = []
    for child in direct:
        # Re-parent structural features onto the synthetic transcript id.
        new_attrs = _rewrite_attributes(
            child.attrs,
            set_parents=[source_gene_id],
        )
        child_lines.append(
            _format_line(
                child,
                seqid=child.seqid,
                start=child.start,
                end=child.end,
                attrs_raw=new_attrs,
            )
        )
    return TranscriptModel(
        transcript_id=source_gene_id,
        seqid=seqid,
        start=start,
        end=end,
        strand=(strand or "").strip(),
        synthetic=True,
        has_cds=_has_cds(direct),
        has_structure=bool(direct),
        lines=[mrna_line, *child_lines],
    )


def dedupe_gff_lines(lines: list[str]) -> list[str]:
    """
    Keep first occurrence of each feature ``ID``; lines without an ID are
    keyed by their full content so identical rows collapse once.
    """
    out: list[str] = []
    seen: set[str] = set()
    for line in lines:
        feat = parse_gff_line(line)
        if feat is None:
            key = line.rstrip("\n")
        elif feat.feature_id:
            key = f"id:{feat.feature_id}"
        else:
            key = f"line:{feat.raw_line}"
        if key in seen:
            continue
        seen.add(key)
        out.append(line if isinstance(line, str) else feat.raw_line)
    return out


def shift_gff_lines(
    lines: list[str],
    *,
    new_seqid: str,
    origin: int,
) -> list[str]:
    """
    Rewrite GFF lines onto ``new_seqid`` with coordinates shifted so
    genomic ``origin`` becomes position 1 (``origin`` is 1-based).

    Fail-closed: non-comment, non-empty lines that do not parse raise
    ``ValueError`` (never silently dropped).
    """
    if origin < 1:
        raise ValueError("origin must be >= 1")
    delta = origin - 1
    out: list[str] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        feat = parse_gff_line(line)
        if feat is None:
            preview = stripped[:120]
            raise ValueError(
                f"unparseable GFF line at index {index}: {preview!r}"
            )
        start = feat.start - delta
        end = feat.end - delta
        if start < 1 or end < start:
            raise ValueError(
                f"feature {feat.feature_id or feat.type} falls outside window"
            )
        out.append(
            _format_line(
                feat,
                seqid=new_seqid,
                start=start,
                end=end,
            )
        )
    return out


def write_fasta(path, *, seqid: str, sequence: str, width: int = 70) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{seqid}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i : i + width])
            fh.write("\n")


def write_gff(path, lines: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for line in lines:
            fh.write(line)
            if not line.endswith("\n"):
                fh.write("\n")
