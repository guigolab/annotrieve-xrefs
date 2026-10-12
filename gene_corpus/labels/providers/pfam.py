"""Pfam label provider (Pfam-A.clans.tsv.gz)."""
from __future__ import annotations

import gzip
import os
import re
import shutil
import ssl
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path

from gene_corpus.labels.providers.base import LabelTerm

PFAM_CLANS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
    "Pfam-A.clans.tsv.gz"
)
CLANS_CACHE_NAME = "Pfam-A.clans.tsv.gz"
MAX_ALIASES_PER_TERM = 5

_PF_ACC_RE = re.compile(r"^(PF\d+)(?:\.\d+)?$", re.IGNORECASE)

_UA = "annotrieve-gene-corpus-labels/1.0 (+https://genome.crg.es/annotrieve)"
_SSL_CTX = ssl._create_unverified_context()


def normalize_pfam_accession(raw: str) -> str | None:
    """Return unversioned PF###### or None if not a Pfam accession."""
    m = _PF_ACC_RE.match(raw.strip())
    if m is None:
        return None
    return m.group(1).upper()


def _urlopen(url: str, *, timeout: float = 600):
    req = urllib.request.Request(
        url, headers={"Accept": "*/*", "User-Agent": _UA}
    )
    return urllib.request.urlopen(req, context=_SSL_CTX, timeout=timeout)


def _download_to(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=dest.name + ".", suffix=".part", dir=str(dest.parent)
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with _urlopen(url) as resp, tmp_path.open("wb") as out:
            shutil.copyfileobj(resp, out, length=1024 * 1024)
        os.replace(tmp_path, dest)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _cap_aliases(aliases: list[str], *, label: str) -> tuple[str, ...]:
    out: list[str] = []
    seen: set[str] = {label}
    for alias in aliases:
        a = alias.strip()
        if not a or a in seen:
            continue
        seen.add(a)
        out.append(a)
        if len(out) >= MAX_ALIASES_PER_TERM:
            break
    return tuple(out)


class PfamLabelProvider:
    """Download and stream-parse Pfam-A clan TSV labels."""

    namespace = "pfam"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache = cache_dir / "pfam"
        cache.mkdir(parents=True, exist_ok=True)
        dest = cache / CLANS_CACHE_NAME

        if not force and dest.is_file() and dest.stat().st_size > 0:
            print(
                f"labels pfam: using cached {dest}",
                file=sys.stderr,
                flush=True,
            )
            return dest

        print(
            f"labels pfam: downloading {PFAM_CLANS_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_to(PFAM_CLANS_URL, dest)
        return dest

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from a .tsv.gz or plain TSV fixture file."""
        if path.is_dir():
            path = path / CLANS_CACHE_NAME
        if not path.is_file():
            return
        yield from _iter_clans_tsv(path)


def _iter_clans_tsv(path: Path) -> Iterator[LabelTerm]:
    # Columns: pfam_acc, clan_acc, clan_id, pfam_id, description
    # clan_* intentionally ignored (shared across families → FTS noise).
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            acc = normalize_pfam_accession(parts[0])
            short_id = parts[3].strip()
            label = parts[4].strip()
            if acc is None or not label:
                continue
            aliases: list[str] = []
            if short_id:
                aliases.append(short_id)
            yield LabelTerm(
                accession=acc,
                label=label,
                aliases=_cap_aliases(aliases, label=label),
            )
