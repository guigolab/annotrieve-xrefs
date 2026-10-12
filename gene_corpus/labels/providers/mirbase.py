"""miRBase label provider (official dumps, RNAcentral fallback)."""
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

MIRBASE_BASES = (
    "https://mirbase.org/download/CURRENT/database_files",
    "https://www.mirbase.org/download/CURRENT/database_files",
)
MIRBASE_OFFICIAL_FILES = (
    "mirna.txt",
    "mirna_mature.txt",
    "dead_mirna.txt",
    "aliases.txt",
)
RNACENTRAL_MIRBASE_TSV_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/"
    "id_mapping/database_mappings/mirbase.tsv"
)
RNACENTRAL_CACHE_NAME = "rnacentral_mirbase.tsv"
SOURCE_MARKER = "source.txt"
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


def _write_source(cache: Path, source: str) -> None:
    (cache / SOURCE_MARKER).write_text(source + "\n", encoding="utf-8")


def _read_source(cache: Path) -> str | None:
    marker = cache / SOURCE_MARKER
    if not marker.is_file():
        return None
    return marker.read_text(encoding="utf-8").strip() or None


def _official_complete(cache: Path) -> bool:
    return all((cache / name).is_file() for name in MIRBASE_OFFICIAL_FILES[:3])


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


class MirbaseLabelProvider:
    """Download and stream-parse miRBase accession labels.

    Prefers official ``database_files`` (mirna / mature / dead / aliases).
    If those hosts are unreachable, falls back to RNAcentral's mirbase.tsv
    id-mapping (accession → name; dead entries may be incomplete).
    """

    namespace = "mirbase"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache = cache_dir / "mirbase"
        cache.mkdir(parents=True, exist_ok=True)

        if not force:
            if _official_complete(cache):
                print(
                    f"labels mirbase: using cached official dumps in {cache}",
                    file=sys.stderr,
                    flush=True,
                )
                _write_source(cache, "official")
                return cache
            rnac = cache / RNACENTRAL_CACHE_NAME
            if rnac.is_file() and rnac.stat().st_size > 0:
                print(
                    f"labels mirbase: using cached RNAcentral fallback {rnac}",
                    file=sys.stderr,
                    flush=True,
                )
                _write_source(cache, "rnacentral")
                return cache

        try:
            self._download_official(cache)
            _write_source(cache, "official")
            return cache
        except Exception as exc:
            print(
                f"labels mirbase: official download failed ({exc}); "
                f"trying RNAcentral fallback",
                file=sys.stderr,
                flush=True,
            )
            self._download_rnacentral(cache)
            _write_source(cache, "rnacentral")
            return cache

    def _download_official(self, cache: Path) -> None:
        last_err: Exception | None = None
        for base in MIRBASE_BASES:
            try:
                for name in MIRBASE_OFFICIAL_FILES:
                    url = f"{base}/{name}"
                    dest = cache / name
                    print(
                        f"labels mirbase: downloading {url}",
                        file=sys.stderr,
                        flush=True,
                    )
                    _download_to(url, dest)
                if not _official_complete(cache):
                    raise RuntimeError("official mirbase files incomplete")
                return
            except Exception as exc:  # noqa: BLE001 — try next mirror
                last_err = exc
                print(
                    f"labels mirbase: mirror failed ({base}): {exc}",
                    file=sys.stderr,
                    flush=True,
                )
        raise RuntimeError(
            f"all mirbase.org mirrors failed: {last_err}"
        ) from last_err

    def _download_rnacentral(self, cache: Path) -> None:
        dest = cache / RNACENTRAL_CACHE_NAME
        print(
            f"labels mirbase: downloading {RNACENTRAL_MIRBASE_TSV_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_to(RNACENTRAL_MIRBASE_TSV_URL, dest)
        if dest.stat().st_size == 0:
            raise RuntimeError("RNAcentral mirbase.tsv is empty")

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from a cache dir, official file set, or RNAcentral TSV."""
        if path.is_file():
            if path.name == RNACENTRAL_CACHE_NAME or path.suffix == ".tsv":
                yield from _iter_rnacentral_tsv(path)
                return
            raise ValueError(f"unsupported mirbase source file: {path}")

        source = _read_source(path)
        rnac = path / RNACENTRAL_CACHE_NAME
        if source == "rnacentral" or (
            source is None and rnac.is_file() and not _official_complete(path)
        ):
            yield from _iter_rnacentral_tsv(rnac)
            return
        yield from _iter_official(path)


def _load_aliases(path: Path) -> dict[str, list[str]]:
    """aliases.txt: accession \\t name1;name2;…"""
    out: dict[str, list[str]] = {}
    if not path.is_file():
        return out
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            acc, rest = line.split("\t", 1)
            acc = acc.strip()
            if not acc.startswith("MI"):
                continue
            names = [p.strip() for p in rest.split(";") if p.strip()]
            if names:
                out[acc] = names
    return out


def _iter_official(cache: Path) -> Iterator[LabelTerm]:
    alias_map = _load_aliases(cache / "aliases.txt")
    seen: set[str] = set()

    mirna = cache / "mirna.txt"
    if mirna.is_file():
        with mirna.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                # auto_mirna, mirna_acc, mirna_id, …
                if len(parts) < 3:
                    continue
                acc = parts[1].strip()
                label = parts[2].strip()
                if not acc.startswith("MI") or acc.startswith("MIMAT"):
                    continue
                if not label or acc in seen:
                    continue
                seen.add(acc)
                yield LabelTerm(
                    accession=acc,
                    label=label,
                    aliases=_cap_aliases(alias_map.get(acc, []), label=label),
                )

    mature = cache / "mirna_mature.txt"
    if mature.is_file():
        with mature.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                # auto_mature, mature_name, mature_acc, …
                if len(parts) < 3:
                    continue
                label = parts[1].strip()
                acc = parts[2].strip()
                if not acc.startswith("MIMAT") or not label or acc in seen:
                    continue
                seen.add(acc)
                yield LabelTerm(
                    accession=acc,
                    label=label,
                    aliases=_cap_aliases(alias_map.get(acc, []), label=label),
                )

    dead = cache / "dead_mirna.txt"
    if dead.is_file():
        with dead.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                acc = parts[0].strip()
                label = parts[1].strip()
                if not acc.startswith("MI") or not label or acc in seen:
                    continue
                seen.add(acc)
                yield LabelTerm(
                    accession=acc,
                    label=label,
                    aliases=_cap_aliases(alias_map.get(acc, []), label=label),
                )


def _iter_rnacentral_tsv(path: Path) -> Iterator[LabelTerm]:
    """RNAcentral mapping: URS, MIRBASE, accession, taxid, type, name."""
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            acc = parts[2].strip()
            label = parts[5].strip()
            if not acc.startswith("MI") or not label or acc in seen:
                continue
            seen.add(acc)
            yield LabelTerm(accession=acc, label=label, aliases=())
