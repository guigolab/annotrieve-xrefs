"""UniProtKB/Swiss-Prot label provider (uniprot_sprot.fasta.gz)."""
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

SWISSPROT_FASTA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot_sprot.fasta.gz"
)
FASTA_CACHE_NAME = "uniprot_sprot.fasta.gz"
MAX_ALIASES_PER_TERM = 5

# >sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 ...
_HEADER_RE = re.compile(
    r"^>(?P<db>sp|tr)\|(?P<ac>[^|]+)\|(?P<entry>\S+)\s+"
    r"(?P<name>.*?)"
    r"(?:\s+OS=|\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)"
)
_GN_RE = re.compile(r"\bGN=(\S+)")

_UA = "annotrieve-gene-corpus-labels/1.0 (+https://genome.crg.es/annotrieve)"
_SSL_CTX = ssl._create_unverified_context()


def _urlopen(url: str, *, timeout: float = 1800):
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


def parse_swissprot_header(line: str) -> LabelTerm | None:
    """Parse one FASTA header line into a LabelTerm, or None if invalid."""
    line = line.rstrip("\n")
    if not line.startswith(">"):
        return None
    m = _HEADER_RE.match(line)
    if m is None:
        # Fallback: minimal >sp|AC|ENTRY rest
        parts = line[1:].split("|", 2)
        if len(parts) < 3 or parts[0] not in {"sp", "tr"}:
            return None
        acc = parts[1].strip()
        rest = parts[2].strip()
        entry, _, after = rest.partition(" ")
        name = after
        for marker in (" OS=", " OX=", " GN=", " PE=", " SV="):
            idx = name.find(marker)
            if idx != -1:
                name = name[:idx]
        name = name.strip()
        if not acc or not name:
            return None
        aliases = [entry] if entry else []
        gn = _GN_RE.search(line)
        if gn:
            aliases.extend(gn.group(1).split())
        return LabelTerm(
            accession=acc,
            label=name,
            aliases=_cap_aliases(aliases, label=name),
        )

    acc = m.group("ac").strip()
    entry = m.group("entry").strip()
    name = m.group("name").strip()
    if not acc or not name:
        return None
    aliases: list[str] = []
    if entry:
        aliases.append(entry)
    gn = _GN_RE.search(line)
    if gn:
        # GN=TP53 or GN=TP53;... — UniProt uses spaces for synonyms after first
        # but typically one token until next field; split on spaces/semicolons.
        for tok in re.split(r"[\s;]+", gn.group(1)):
            if tok:
                aliases.append(tok)
    return LabelTerm(
        accession=acc,
        label=name,
        aliases=_cap_aliases(aliases, label=name),
    )


class SwissprotLabelProvider:
    """Stream UniProt Swiss-Prot FASTA headers into label terms."""

    namespace = "uniprot/swiss-prot"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache = cache_dir / "uniprot_swiss_prot"
        cache.mkdir(parents=True, exist_ok=True)
        dest = cache / FASTA_CACHE_NAME

        if not force and dest.is_file() and dest.stat().st_size > 0:
            print(
                f"labels uniprot/swiss-prot: using cached {dest}",
                file=sys.stderr,
                flush=True,
            )
            return dest

        print(
            f"labels uniprot/swiss-prot: downloading {SWISSPROT_FASTA_URL}",
            file=sys.stderr,
            flush=True,
        )
        _download_to(SWISSPROT_FASTA_URL, dest)
        return dest

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield terms from a .fasta.gz or plain .fasta / fixture file."""
        if path.is_dir():
            path = path / FASTA_CACHE_NAME
        if not path.is_file():
            return
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
            for line in fh:
                if not line.startswith(">"):
                    continue
                term = parse_swissprot_header(line)
                if term is not None:
                    yield term
