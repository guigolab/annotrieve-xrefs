"""Resolve on-disk annotation GFF paths and map GFF seqids to INSDC accessions."""
from __future__ import annotations

import json
import re
import threading
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote

from settings import annotations_dir, settings

CHROMOSOMES_FILENAME = "chromosomes.json"
CHR_ALIASES_FILENAME = "chr_aliases.tsv"

_ASSEMBLY_RE = re.compile(r"^(GC[AF])_(\d+)(?:\.\d+)?$", re.IGNORECASE)


class _MappedWithoutAccession(Exception):
    """chromosomes/aliases matched a row that has empty genbank+refseq fields."""

_CACHE_LOCK = threading.Lock()


@dataclass(frozen=True)
class AnnotationFiles:
    """Resolved paths for one annotation's on-disk artifacts."""

    gff_path: Path
    csi_path: Path
    taxid: int
    assembly_accession: str
    assembly_dir: Path
    source_database: str
    # Other folders that may hold chromosomes.json / chr_aliases.tsv: Annotrieve leaves them in the
    # annotation's current taxid folder while the GFF is still under a retired taxid folder.
    extra_dirs: tuple[Path, ...] = ()


class AnnotationFilesError(Exception):
    """Local file / seqid resolution failure (maps to API error codes)."""

    def __init__(self, code: str, message: str, **extra) -> None:
        super().__init__(message)
        self.code = code
        self.extra = extra


def annotations_root() -> Path:
    root = annotations_dir()
    if root is None:
        raise AnnotationFilesError(
            "annotations_unavailable",
            "ANNOTATIONS_DIR is not configured",
        )
    if not root.is_dir():
        raise AnnotationFilesError(
            "annotations_unavailable",
            f"annotations directory not found: {root}",
        )
    return root


def _rebase_gff_path(gff_path: str, root: Path) -> Path | None:
    """Map harvest-time absolute gff_path onto ANNOTATIONS_DIR when possible."""
    raw = (gff_path or "").strip()
    if not raw:
        return None
    path = Path(raw)
    prefix = (settings.ANNOTATIONS_GFF_PATH_PREFIX or "").strip()
    if prefix:
        pref = Path(prefix)
        try:
            relative = path.relative_to(pref)
            candidate = root / relative
            if candidate.is_file():
                return candidate
        except ValueError:
            pass
        # Also accept prefix as a string prefix (symlink / trailing slash variants).
        pref_s = str(pref).rstrip("/") + "/"
        raw_s = raw if raw.startswith("/") else f"/{raw}"
        if raw_s.startswith(pref_s) or raw.startswith(str(pref).rstrip("/") + "/"):
            tail = raw[len(str(pref).rstrip("/")) :].lstrip("/")
            candidate = root / tail
            if candidate.is_file():
                return candidate

    # If gff_path already lives under root, use it.
    try:
        resolved = path.resolve()
        root_res = root.resolve()
        if root_res in resolved.parents or resolved == root_res:
            if resolved.is_file():
                return resolved
    except OSError:
        pass

    # Strip a leading absolute root-like segment and try the basename path
    # under taxid/assembly via the convention caller.
    return None


def _tail_gff_path(root: Path, gff_path: str) -> Path | None:
    """``root/<folder>/<assembly>/<file>`` from the last three parts of the stored path.

    The folder on disk is not always the annotation's current taxid (NCBI re-assigns taxids and
    Annotrieve used to leave the files in the old folder), but the stored path carries the real one.
    """
    parts = Path((gff_path or "").strip()).parts[-3:]
    if len(parts) < 3 or any(p in ("", "/", "..", ".") for p in parts):
        return None
    candidate = root.joinpath(*parts)
    return candidate if candidate.is_file() else None


def _convention_gff_path(
    root: Path,
    *,
    taxid: int,
    assembly_accession: str,
    gff_path: str,
) -> Path | None:
    name = Path((gff_path or "").strip()).name
    if not name:
        return None
    candidate = root / str(taxid) / assembly_accession / name
    if candidate.is_file():
        return candidate
    # Some community uploads sit directly under the taxid dir.
    alt = root / str(taxid) / name
    if alt.is_file():
        return alt
    return None


def resolve_annotation_files(
    *,
    taxid: int,
    assembly_accession: str,
    gff_path: str,
    source_database: str = "",
) -> AnnotationFiles:
    """
    Locate the bgzipped GFF + CSI for one annotation.

    Resolution order:
      1. Rebase ``gff_path`` onto ``ANNOTATIONS_DIR`` (optional harvest prefix).
      2. Convention ``ANNOTATIONS_DIR/{taxid}/{assembly}/{basename(gff_path)}``.
    """
    root = annotations_root()
    resolved = _rebase_gff_path(gff_path, root)
    if resolved is None:
        resolved = _tail_gff_path(root, gff_path)
    if resolved is None:
        resolved = _convention_gff_path(
            root,
            taxid=taxid,
            assembly_accession=assembly_accession,
            gff_path=gff_path,
        )
    if resolved is None or not resolved.is_file():
        raise AnnotationFilesError(
            "annotation_file_not_found",
            "annotation GFF file not found on disk",
            taxid=taxid,
            assembly_accession=assembly_accession,
            gff_basename=Path((gff_path or "").strip()).name or None,
        )
    csi = Path(str(resolved) + ".csi")
    if not csi.is_file():
        # pysam also accepts .tbi
        tbi = Path(str(resolved) + ".tbi")
        if tbi.is_file():
            csi = tbi
        else:
            raise AnnotationFilesError(
                "annotation_file_not_found",
                "annotation GFF index (.csi/.tbi) not found on disk",
                taxid=taxid,
                assembly_accession=assembly_accession,
                gff_path=str(resolved),
            )
    own = root / str(int(taxid)) / str(assembly_accession)
    extra = (own,) if own != resolved.parent and own.is_dir() else ()
    return AnnotationFiles(
        gff_path=resolved,
        csi_path=csi,
        taxid=int(taxid),
        assembly_accession=str(assembly_accession),
        assembly_dir=resolved.parent,
        source_database=source_database or "",
        extra_dirs=extra,
    )


def _sibling_assembly_dirs(taxid_dir: Path, assembly_accession: str) -> list[Path]:
    """Paired GCA/GCF dirs under the same taxid that share the numeric id."""
    m = _ASSEMBLY_RE.match((assembly_accession or "").strip())
    if not m or not taxid_dir.is_dir():
        return []
    digits = m.group(2)
    other = "GCF" if m.group(1).upper() == "GCA" else "GCA"
    out: list[Path] = []
    try:
        for child in taxid_dir.iterdir():
            if not child.is_dir():
                continue
            cm = _ASSEMBLY_RE.match(child.name)
            if cm and cm.group(2) == digits and cm.group(1).upper() == other:
                out.append(child)
    except OSError:
        return []
    return sorted(out)


def _candidates(files: AnnotationFiles, filename: str) -> list[Path]:
    cands: list[Path] = []
    for ass_dir in (files.assembly_dir, *files.extra_dirs):
        cands.append(ass_dir / filename)
        for sibling in _sibling_assembly_dirs(ass_dir.parent, files.assembly_accession):
            cands.append(sibling / filename)
    return cands


def _chromosomes_candidates(files: AnnotationFiles) -> list[Path]:
    return _candidates(files, CHROMOSOMES_FILENAME)


def _aliases_candidates(files: AnnotationFiles) -> list[Path]:
    return _candidates(files, CHR_ALIASES_FILENAME)


@lru_cache(maxsize=256)
def _load_chromosomes(path_str: str) -> tuple[dict, ...]:
    path = Path(path_str)
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, list):
        return ()
    return tuple(row for row in data if isinstance(row, dict))


@lru_cache(maxsize=256)
def _load_aliases_index(path_str: str) -> dict[str, str]:
    """alias → refName (first column of chr_aliases.tsv)."""
    index: dict[str, str] = {}
    with open(path_str, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            ref = parts[0].strip()
            if not ref:
                continue
            for alias in parts:
                alias = alias.strip()
                if alias:
                    index[alias] = ref
    return index


def _chromosome_identifiers(row: dict) -> set[str]:
    ids: set[str] = set()
    for key in (
        "chr_name",
        "sequence_name",
        "ucsc_style_name",
        "genbank_accession",
        "refseq_accession",
        "canonical_id",
    ):
        val = row.get(key)
        if not val:
            continue
        v = str(val).strip()
        if not v:
            continue
        ids.add(v)
        if "|" in v:
            ids.update(part.strip() for part in v.split("|") if part.strip())
        if "." in v:
            ids.add(v.split(".", 1)[0])
    return ids


def _prefer_refseq(source_database: str) -> bool:
    db = (source_database or "").strip().casefold()
    return db in {"refseq", "ncbi"}


def _pick_accession(row: dict, *, prefer_refseq: bool) -> str | None:
    """Return preferred non-empty genbank/refseq field (no format filter)."""
    rs = str(row.get("refseq_accession") or "").strip()
    gb = str(row.get("genbank_accession") or "").strip()
    order = (rs, gb) if prefer_refseq else (gb, rs)
    for cand in order:
        if cand:
            return cand
    return None


@lru_cache(maxsize=4096)
def _resolve_seqid_accession_cached(
    assembly_dir: str,
    source_database: str,
    seqid: str,
    extra_dirs: tuple[str, ...] = (),
) -> str:
    """
    Cached resolver core.

    Returns a mapped accession, an alias ``refName``, or the raw ``seqid``
    (pass-through). Raises ``_MappedWithoutAccession`` when a chromosome row
    matches but both accession fields are empty.
    """
    ass_dir = Path(assembly_dir)
    files = AnnotationFiles(
        gff_path=ass_dir / "_unused.gff.gz",
        csi_path=ass_dir / "_unused.gff.gz.csi",
        taxid=0,
        assembly_accession=ass_dir.name,
        assembly_dir=ass_dir,
        source_database=source_database,
        extra_dirs=tuple(Path(d) for d in extra_dirs),
    )
    prefer_refseq = _prefer_refseq(source_database)

    # Match via chromosomes.json
    for cand in _chromosomes_candidates(files):
        if not cand.is_file():
            continue
        try:
            rows = _load_chromosomes(str(cand.resolve()))
        except (OSError, json.JSONDecodeError, ValueError):
            continue
        for row in rows:
            if seqid in _chromosome_identifiers(row):
                acc = _pick_accession(row, prefer_refseq=prefer_refseq)
                if acc:
                    return acc
                raise _MappedWithoutAccession(seqid)

    # Match via chr_aliases.tsv → refName → chromosomes.json
    ref_name: str | None = None
    for cand in _aliases_candidates(files):
        if not cand.is_file():
            continue
        try:
            index = _load_aliases_index(str(cand.resolve()))
        except OSError:
            continue
        ref_name = index.get(seqid)
        if ref_name:
            break
    if ref_name:
        for cand in _chromosomes_candidates(files):
            if not cand.is_file():
                continue
            try:
                rows = _load_chromosomes(str(cand.resolve()))
            except (OSError, json.JSONDecodeError, ValueError):
                continue
            for row in rows:
                identifiers = _chromosome_identifiers(row)
                if ref_name in identifiers or seqid in identifiers:
                    acc = _pick_accession(row, prefer_refseq=prefer_refseq)
                    if acc:
                        return acc
                    raise _MappedWithoutAccession(seqid)
        return ref_name

    # No chromosomes/aliases hit — pass the GFF seqid through to NCBI.
    return seqid


def resolve_seqid_accession(
    files: AnnotationFiles,
    seqid: str,
) -> str:
    """
    Map a GFF seqid to an accession for NCBI efetch.

    Looks up ``chromosomes.json`` / ``chr_aliases.tsv``. Unmapped seqids are
    returned as-is (pass-through). Raises ``AnnotationFilesError`` only for
    empty seqid or a matched chromosome row with empty genbank+refseq fields.
    """
    raw = unquote((seqid or "").strip())
    if not raw:
        raise AnnotationFilesError(
            "seqid_not_resolved",
            "empty seqid",
            seqid=seqid,
        )

    try:
        accession = _resolve_seqid_accession_cached(
            str(files.assembly_dir.resolve()),
            files.source_database or "",
            raw,
            tuple(str(d.resolve()) for d in files.extra_dirs),
        )
        if accession == raw and raw.isdigit():
            # Unmapped bare number (Ensembl "1", "17"): NCBI would read it as a GI record id and
            # could return an unrelated sequence, so it must not be passed through.
            raise AnnotationFilesError(
                "seqid_not_resolved",
                f"seqid {raw!r} is a bare number with no chromosome mapping for this assembly "
                "(NCBI would treat it as a GI id)",
                seqid=raw,
                assembly_accession=files.assembly_accession,
                taxid=files.taxid,
            )
        return accession
    except _MappedWithoutAccession:
        raise AnnotationFilesError(
            "seqid_not_resolved",
            (
                f"could not resolve seqid {raw!r} to an accession "
                "(chromosome row has no genbank/refseq)"
            ),
            seqid=raw,
            assembly_accession=files.assembly_accession,
            taxid=files.taxid,
        ) from None
    except OSError as exc:
        raise AnnotationFilesError(
            "seqid_not_resolved",
            f"could not resolve seqid {raw!r}: {exc}",
            seqid=raw,
            assembly_accession=files.assembly_accession,
            taxid=files.taxid,
        ) from exc


def clear_annotation_file_caches() -> None:
    """Drop LRU caches (tests)."""
    with _CACHE_LOCK:
        _load_chromosomes.cache_clear()
        _load_aliases_index.cache_clear()
        _resolve_seqid_accession_cached.cache_clear()
