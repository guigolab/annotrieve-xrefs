"""Rfam label provider (family + dead_family MySQL dumps)."""
from __future__ import annotations

import gzip
import os
import shutil
import ssl
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path

from gene_corpus.labels.providers.base import LabelTerm

RFAM_FAMILY_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/"
    "database_files/family.txt.gz"
)
RFAM_DEAD_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/"
    "database_files/dead_family.txt.gz"
)
FAMILY_NAME = "family.txt"
DEAD_NAME = "dead_family.txt"
MAX_ALIASES_PER_TERM = 5

_UA = "annotrieve-gene-corpus-labels/1.0 (+https://genome.crg.es/annotrieve)"
_SSL_CTX = ssl._create_unverified_context()


def _urlopen(url: str, *, timeout: float = 600):
    req = urllib.request.Request(
        url, headers={"Accept": "*/*", "User-Agent": _UA}
    )
    return urllib.request.urlopen(req, context=_SSL_CTX, timeout=timeout)


def _download_gunzip(url: str, dest_txt: Path) -> None:
    """Download a .gz URL and write the decompressed text file atomically."""
    dest_txt.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=dest_txt.name + ".", suffix=".part", dir=str(dest_txt.parent)
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with _urlopen(url) as resp, gzip.GzipFile(fileobj=resp) as gz, tmp_path.open(
            "wb"
        ) as out:
            shutil.copyfileobj(gz, out, length=1024 * 1024)
        os.replace(tmp_path, dest_txt)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _cap_aliases(aliases: list[str], *, label: str) -> tuple[str, ...]:
    out: list[str] = []
    seen: set[str] = {label}
    for alias in aliases:
        a = alias.strip()
        if not a or a.upper() == "NULL" or a in seen:
            continue
        seen.add(a)
        out.append(a)
        if len(out) >= MAX_ALIASES_PER_TERM:
            break
    return tuple(out)


def _cell(parts: list[str], idx: int) -> str:
    if idx >= len(parts):
        return ""
    return parts[idx].strip()


def _previous_id_tokens(raw: str) -> list[str]:
    if not raw or raw.upper() == "NULL":
        return []
    # previous_id is a tinytext of historical short names (often ';'-separated).
    return [p.strip() for p in raw.replace(",", ";").split(";") if p.strip()]


class RfamLabelProvider:
    """Download and stream-parse Rfam family labels (incl. dead families)."""

    namespace = "rfam"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache = cache_dir / "rfam"
        cache.mkdir(parents=True, exist_ok=True)
        family = cache / FAMILY_NAME
        dead = cache / DEAD_NAME

        if (
            not force
            and family.is_file()
            and family.stat().st_size > 0
            and dead.is_file()
            and dead.stat().st_size > 0
        ):
            print(
                f"labels rfam: using cached dumps in {cache}",
                file=sys.stderr,
                flush=True,
            )
            return cache

        print(
            f"labels rfam: downloading {RFAM_FAMILY_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_gunzip(RFAM_FAMILY_URL, family)
        print(
            f"labels rfam: downloading {RFAM_DEAD_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_gunzip(RFAM_DEAD_URL, dead)
        return cache

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from a cache dir or a single family TSV fixture file."""
        if path.is_file():
            yield from _iter_family_file(path)
            return
        family = path / FAMILY_NAME
        dead = path / DEAD_NAME
        seen: set[str] = set()
        if family.is_file():
            for term in _iter_family_file(family):
                seen.add(term.accession)
                yield term
        if dead.is_file():
            for term in _iter_dead_file(dead):
                if term.accession in seen:
                    continue
                seen.add(term.accession)
                yield term


def _iter_family_file(path: Path) -> Iterator[LabelTerm]:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            acc = _cell(parts, 0)
            label = _cell(parts, 1)
            if not acc.startswith("RF") or not label:
                continue
            desc = _cell(parts, 3)
            prev = _cell(parts, 10)
            alias_bits: list[str] = []
            if desc and desc != label:
                alias_bits.append(desc)
            alias_bits.extend(_previous_id_tokens(prev))
            yield LabelTerm(
                accession=acc,
                label=label,
                aliases=_cap_aliases(alias_bits, label=label),
            )


def _iter_dead_file(path: Path) -> Iterator[LabelTerm]:
    # dead_family: rfam_acc, rfam_id, comment, forward_to, title, user
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            acc = _cell(parts, 0)
            label = _cell(parts, 1)
            if not acc.startswith("RF") or not label:
                continue
            title = _cell(parts, 4)
            # forward_to (col 3) intentionally ignored — not an alt_id.
            yield LabelTerm(
                accession=acc,
                label=label,
                aliases=_cap_aliases([title] if title else [], label=label),
            )
