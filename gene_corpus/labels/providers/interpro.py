"""InterPro label provider (entry.list + short_names.dat)."""
from __future__ import annotations

import os
import shutil
import ssl
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path

from gene_corpus.labels.providers.base import LabelTerm

ENTRY_LIST_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
)
SHORT_NAMES_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/short_names.dat"
)
ENTRY_LIST_NAME = "entry.list"
SHORT_NAMES_NAME = "short_names.dat"
MAX_ALIASES_PER_TERM = 5

_UA = "annotrieve-gene-corpus-labels/1.0 (+https://genome.crg.es/annotrieve)"
_SSL_CTX = ssl._create_unverified_context()


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


class InterproLabelProvider:
    """Download and stream-parse InterPro entry labels."""

    namespace = "interpro"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache = cache_dir / "interpro"
        cache.mkdir(parents=True, exist_ok=True)
        entry = cache / ENTRY_LIST_NAME
        shorts = cache / SHORT_NAMES_NAME

        if (
            not force
            and entry.is_file()
            and entry.stat().st_size > 0
            and shorts.is_file()
            and shorts.stat().st_size > 0
        ):
            print(
                f"labels interpro: using cached dumps in {cache}",
                file=sys.stderr,
                flush=True,
            )
            return cache

        print(
            f"labels interpro: downloading {ENTRY_LIST_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_to(ENTRY_LIST_URL, entry)
        print(
            f"labels interpro: downloading {SHORT_NAMES_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_to(SHORT_NAMES_URL, shorts)
        return cache

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from a cache dir (or entry.list file for tests)."""
        if path.is_file():
            yield from _iter_entry_list(path, short_names={})
            return
        entry = path / ENTRY_LIST_NAME
        shorts_path = path / SHORT_NAMES_NAME
        short_names = _load_short_names(shorts_path) if shorts_path.is_file() else {}
        if entry.is_file():
            yield from _iter_entry_list(entry, short_names=short_names)


def _load_short_names(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            acc, short = line.split("\t", 1)
            acc = acc.strip()
            short = short.strip()
            if acc.startswith("IPR") and short:
                out[acc] = short
    return out


def _iter_entry_list(
    path: Path, *, short_names: dict[str, str]
) -> Iterator[LabelTerm]:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if i == 0 and parts[0].strip().upper() == "ENTRY_AC":
                continue
            # ENTRY_AC, ENTRY_TYPE, ENTRY_NAME — type intentionally ignored.
            acc = parts[0].strip() if parts else ""
            label = parts[2].strip() if len(parts) > 2 else ""
            if not acc.startswith("IPR") or not label:
                continue
            alias_bits: list[str] = []
            short = short_names.get(acc)
            if short:
                alias_bits.append(short)
            yield LabelTerm(
                accession=acc,
                label=label,
                aliases=_cap_aliases(alias_bits, label=label),
            )
