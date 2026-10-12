"""Gene Ontology label provider (go-basic.obo)."""
from __future__ import annotations

import os
import re
import shutil
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path

from gene_corpus.labels.providers.base import LabelTerm
from gene_corpus.stream.emit_util import normalize_go_token

GO_BASIC_OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GO_CACHE_NAME = "go-basic.obo"

# Cap aliases stored per term; prefer EXACT then RELATED when truncating.
MAX_ALIASES_PER_TERM = 5
_SYNONYM_SCOPE_RANK = {
    "EXACT": 0,
    "RELATED": 1,
    "NARROW": 2,
    "BROAD": 3,
}
_SYNONYM_RE = re.compile(
    r'^synonym:\s*"((?:\\.|[^"\\])*)"\s*([A-Z_]+)?'
)


def _select_aliases(
    scoped: list[tuple[str, str]],
    *,
    limit: int = MAX_ALIASES_PER_TERM,
) -> list[str]:
    """
    Keep up to *limit* synonym strings, preferring EXACT then RELATED.

    *scoped* is ``(text, scope)`` in OBO order. Stable within each scope.
    """
    if not scoped:
        return []
    ranked = sorted(
        enumerate(scoped),
        key=lambda item: (
            _SYNONYM_SCOPE_RANK.get(item[1][1], 99),
            item[0],
        ),
    )
    out: list[str] = []
    seen: set[str] = set()
    for _, (text, _scope) in ranked:
        if text in seen:
            continue
        seen.add(text)
        out.append(text)
        if len(out) >= limit:
            break
    return out


class GoLabelProvider:
    """Download and stream-parse Gene Ontology basic OBO."""

    namespace = "go"

    def download(self, cache_dir: Path, *, force: bool = False) -> Path:
        cache_dir.mkdir(parents=True, exist_ok=True)
        dest = cache_dir / GO_CACHE_NAME
        if dest.is_file() and not force:
            return dest

        print(
            f"labels go: downloading {GO_BASIC_OBO_URL}",
            file=sys.stderr,
            flush=True,
        )
        fd, tmp_name = tempfile.mkstemp(
            prefix="go-basic.", suffix=".obo", dir=str(cache_dir)
        )
        os.close(fd)
        tmp_path = Path(tmp_name)
        try:
            with urllib.request.urlopen(GO_BASIC_OBO_URL) as resp, tmp_path.open(
                "wb"
            ) as out:
                shutil.copyfileobj(resp, out, length=1024 * 1024)
            os.replace(tmp_path, dest)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise
        return dest

    def iter_terms(self, path: Path) -> Iterator[LabelTerm]:
        """Yield GO terms from an OBO file (stanza stream), including obsolete."""
        term_lines: list[str] = []
        in_term = False
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("["):
                    if in_term and term_lines:
                        term = _parse_term_stanza(term_lines)
                        if term is not None:
                            yield term
                    term_lines = []
                    in_term = line.strip() == "[Term]"
                    continue
                if in_term:
                    if line.strip() == "":
                        term = _parse_term_stanza(term_lines)
                        if term is not None:
                            yield term
                        term_lines = []
                        in_term = False
                    else:
                        term_lines.append(line)
            if in_term and term_lines:
                term = _parse_term_stanza(term_lines)
                if term is not None:
                    yield term


def _parse_term_stanza(lines: list[str]) -> LabelTerm | None:
    term_id: str | None = None
    name: str | None = None
    alt_ids: list[str] = []
    scoped_aliases: list[tuple[str, str]] = []
    seen_alias: set[str] = set()

    for raw in lines:
        if raw.startswith("id:"):
            term_id = normalize_go_token(raw[3:].strip())
        elif raw.startswith("name:"):
            name = raw[5:].strip()
        elif raw.startswith("alt_id:"):
            alt = normalize_go_token(raw[7:].strip())
            if alt is not None:
                alt_ids.append(alt)
        elif raw.startswith("synonym:"):
            m = _SYNONYM_RE.match(raw)
            if m is None:
                continue
            syn = m.group(1).replace('\\"', '"').strip()
            if not syn or syn in seen_alias:
                continue
            seen_alias.add(syn)
            scope = (m.group(2) or "").strip().upper() or "RELATED"
            scoped_aliases.append((syn, scope))
        # is_obsolete is ignored: corpus accessions often still use retired IDs.

    if term_id is None or not name:
        return None
    uniq_alts: list[str] = []
    seen_alt: set[str] = set()
    for alt in alt_ids:
        if alt == term_id or alt in seen_alt:
            continue
        seen_alt.add(alt)
        uniq_alts.append(alt)
    aliases = _select_aliases(scoped_aliases, limit=MAX_ALIASES_PER_TERM)
    return LabelTerm(
        accession=term_id,
        label=name,
        aliases=tuple(aliases),
        alt_ids=tuple(uniq_alts),
    )
