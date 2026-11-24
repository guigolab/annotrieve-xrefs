"""``/hits`` and ``/names`` over presence.sqlite, with alias as a suggestion (API-DESIGN 3)."""
from __future__ import annotations

import json
import threading
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from helpers.cursor import QueryError
from helpers.names import (
    parse_prefix_query,
    parse_taxid_list,
)
from helpers.name_category import classify_name
from helpers.pagination import MAX_NAMES_GET, MAX_PREFIX_EXPAND
from helpers.postings import iter_bits, to_bitset
from helpers.presence import _shared, scope_mask
from helpers.source_codes import SOURCE_NAME_TO_CODE, SOURCES_LEGEND
from services.terms import parse_terms, resolve_terms, term_echo
from settings import settings

ALIAS_MODES = ("off", "suggest", "include")
HITS_FIELDS = [
    "annotation_id", "taxid", "source", "matched", "release_date", "assembly_name",
    "alias_only", "organism_name", "assembly_accession",
]
HITS_SOURCES = {code: {"label": label} for code, label in SOURCES_LEGEND.items()}
_MAX_VARS = 900


def parse_alias_mode(value: str | None) -> str:
    mode = (value or "suggest").strip().lower()
    if mode not in ALIAS_MODES:
        raise QueryError(
            "alias must be one of off, suggest, include",
            code="invalid_alias_mode", value=value,
        )
    return mode


SHAPES = ("compact", "records")
SOURCE_LABELS = {int(code): label for code, label in SOURCES_LEGEND.items()}


def parse_shape(value: str | None) -> str:
    shape = (value or "compact").strip().lower()
    if shape not in SHAPES:
        raise QueryError("shape must be one of compact, records", code="invalid_shape", value=value)
    return shape


def parse_sources(value) -> list[int]:
    """Source codes from a CSV string or a list of names (Ensembl, GenBank, RefSeq); [] = all."""
    if value is None:
        return []
    items = value.split(",") if isinstance(value, str) else list(value)
    codes: list[int] = []
    for item in items:
        if not isinstance(item, str):
            raise QueryError("sources must be names", code="invalid_source", value=str(item))
        key = item.strip().casefold()
        if not key:
            continue
        if key not in SOURCE_NAME_TO_CODE:
            raise QueryError(
                f"unknown source {item.strip()!r}; use Ensembl, GenBank or RefSeq",
                code="invalid_source", value=item.strip(), allowed=sorted(SOURCES_LEGEND.values()))
        code = SOURCE_NAME_TO_CODE[key]
        if code not in codes:
            codes.append(code)
    return codes


def parse_min_matched(value, n_terms: int) -> int | None:
    if value is None or value == "":
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise QueryError("min_matched must be an integer", code="invalid_min_matched",
                         value=str(value)) from None
    if number < 1 or number > n_terms:
        raise QueryError(f"min_matched must be between 1 and the number of terms ({n_terms})",
                         code="invalid_min_matched", value=number, n_terms=n_terms)
    return number


def download_url(gff_path: str | None) -> str | None:
    """Whole-file link on Annotrieve's file server: the last three parts of the stored path."""
    parts = [p for p in str(gff_path or "").split("/") if p][-3:]
    if len(parts) < 3:
        return None
    return f"{settings.ANNOTRIEVE_FILES_BASE.rstrip('/')}/{'/'.join(parts)}"


def _chunks(items: list, size: int = _MAX_VARS) -> Iterable[list]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def _optional_text(raw) -> str | None:
    text = None if raw is None else str(raw).strip()
    return text or None


def _split_by_mode(sym: int, alias: int, mode: str) -> tuple[int, int, int]:
    """Return (matched mask, alias-only mask, alias-extra mask) for one name/term."""
    extra = alias & ~sym
    if mode == "include":
        return sym | alias, 0, extra
    if mode == "off":
        return sym, 0, 0
    return sym, extra, extra


class AnnStore:
    """
    Everything /hits needs to say about each annotation, loaded once per data version.

    ``rec[ord]`` is ``(annotation_id, taxid, source, release_date, assembly_name, organism_name,
    assembly_accession, gff_path)``; ``rank[ord]`` is the position of its id in sorted order (the tie
    breaker of the result order). The JSON text of each annotation is written once, the first time
    it is needed, in pieces around the parts that change per request (the matched term lists), so a
    response is mostly string joining.
    """

    def __init__(self, rows: Iterable) -> None:
        self.rec: dict[int, tuple] = {}
        for r in rows:
            self.rec[int(r["ord"])] = (
                str(r["annotation_id"]), int(r["taxid"]), int(r["source"]),
                _optional_text(r["release_date"]), _optional_text(r["assembly_name"]),
                _optional_text(r["organism_name"]), _optional_text(r["assembly_accession"]),
                r["gff_path"])
        order = sorted(self.rec, key=lambda o: self.rec[o][0])
        self.rank = {o: i for i, o in enumerate(order)}
        self._compact: dict[int, tuple[str, str, str]] = {}
        self._records: dict[int, tuple[str, str]] = {}

    def compact_pieces(self, o: int) -> tuple[str, str, str]:
        pieces = self._compact.get(o)
        if pieces is None:
            aid, taxid, source, date, asm, org, acc, _gff = self.rec[o]
            pieces = (f"[{_j(aid)},{taxid},{source},", f",{_j(date)},{_j(asm)},",
                      f",{_j(org)},{_j(acc)}]")
            self._compact[o] = pieces
        return pieces

    def record_pieces(self, o: int) -> tuple[str, str]:
        pieces = self._records.get(o)
        if pieces is None:
            aid, taxid, source, date, asm, org, acc, gff = self.rec[o]
            pieces = (
                f'{{"annotation_id":{_j(aid)},"organism_name":{_j(org)},"taxid":{taxid},'
                f'"assembly_accession":{_j(acc)},"assembly_name":{_j(asm)},'
                f'"source":{_j(SOURCE_LABELS.get(source))},"release_date":{_j(date)},"n_matched":',
                f',"download_url":{_j(download_url(gff))}}}')
            self._records[o] = pieces
        return pieces


def _j(value) -> str:
    return json.dumps(value, ensure_ascii=False)


_LISTS: dict[tuple, str] = {}


def _ints(values: Sequence[int]) -> str:
    """JSON text of a short list of ints (memoized: most rows share the same few lists)."""
    key = tuple(values)
    text = _LISTS.get(key)
    if text is None:
        text = "[" + ",".join(map(str, key)) + "]"
        if len(_LISTS) < 4096:
            _LISTS[key] = text
    return text


_STORE_LOCK = threading.Lock()


def ann_store(reader) -> AnnStore:
    """The annotation store of the reader's data version (one copy for all threads)."""
    entry = _shared(reader)
    store = entry.get("store")
    if store is None:
        with _STORE_LOCK:
            store = entry.get("store")
            if store is None:
                store = AnnStore(reader.conn.execute(
                    "SELECT ord, annotation_id, taxid, organism_name, assembly_accession, "
                    "assembly_name, source, release_date, gff_path FROM annotation"))
                entry["store"] = store
    return store


@dataclass
class _Search:
    """Result of the term and posting work, before any output format is chosen."""
    base: dict
    n_terms: int
    ann_matched: dict[int, list[int]]
    ann_alias: dict[int, list[int]]
    ords: list[int]          # selected annotations in answer order
    facets: dict
    store: AnnStore


def _search(
    reader,
    names: str | Sequence[str] | None = None,
    *,
    terms: Sequence | None = None,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
    max_names: int = MAX_NAMES_GET,
    sources=None,
    min_matched=None,
) -> _Search:
    mode = parse_alias_mode(alias)
    source_codes = parse_sources(sources)
    raws = parse_terms(names, terms, max_terms=max_names)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")

    conn = reader.conn
    resolved = resolve_terms(conn, raws, mode)
    # Several terms can share a name (``seleno*`` and ``selenow``): index name_id -> terms.
    id_to_terms: dict[int, list[int]] = {}
    terms_out: list[dict] = []
    for t in resolved:
        for nid in t["name_ids"]:
            id_to_terms.setdefault(nid, []).append(t["index"])
        terms_out.append({**term_echo(t), "namespaces": [], "n_hits": 0, "n_alias_extra": 0})
    scope = scope_mask(conn, reader.all_mask, include, exclude)

    term_sym = [0] * len(resolved)
    term_alias = [0] * len(resolved)
    for chunk in _chunks(list(id_to_terms)):
        ph = ",".join("?" for _ in chunk)
        for row in conn.execute(
                f"SELECT name_id, sym, alias FROM posting WHERE name_id IN ({ph})", chunk):
            sym = to_bitset(row["sym"]) & scope
            als = to_bitset(row["alias"]) & scope
            for i in id_to_terms[int(row["name_id"])]:
                term_sym[i] |= sym
                term_alias[i] |= als

    ann_matched: dict[int, list[int]] = {}
    ann_alias: dict[int, list[int]] = {}
    for i, out in enumerate(terms_out):
        sym, als = term_sym[i], term_alias[i]
        matched, alias_only, extra = _split_by_mode(sym, als, mode)
        ns = []
        if sym:
            ns.append("symbol")
        if als and mode != "off":
            ns.append("alias")
        out["namespaces"] = ns
        out["n_hits"] = matched.bit_count()
        out["n_alias_extra"] = extra.bit_count()
        for o in iter_bits(matched):
            ann_matched.setdefault(o, []).append(i)
        for o in iter_bits(alias_only):
            ann_alias.setdefault(o, []).append(i)

    n_terms = len(terms_out)
    needed = parse_min_matched(min_matched, n_terms)
    base = {
        "terms": terms_out, "sources": dict(HITS_SOURCES), "alias": mode,
        "include_taxids": include, "exclude_taxids": exclude,
        "source_filter": [SOURCE_LABELS[c] for c in source_codes], "min_matched": needed,
        "data_version": reader.data_version,
    }
    ords = sorted(set(ann_matched) | set(ann_alias))
    if needed:
        ords = [o for o in ords if len(ann_matched.get(o, ())) >= needed]
    store = ann_store(reader)
    rec = store.rec
    ords = [o for o in ords if o in rec]
    # Counts per source are taken before the sources filter, so a client can show the alternatives.
    facets = {"source": {SOURCE_LABELS[c]: 0 for c in sorted(SOURCE_LABELS)}}
    for o in ords:
        label = SOURCE_LABELS.get(rec[o][2])
        if label:
            facets["source"][label] += 1
    if source_codes:
        ords = [o for o in ords if rec[o][2] in source_codes]
    rank = store.rank
    ords.sort(key=lambda o: (-len(ann_matched.get(o, ())), rank[o]))
    return _Search(base, n_terms, ann_matched, ann_alias, ords, facets, store)


def list_name_hits_presence(
    reader,
    names: str | Sequence[str] | None = None,
    *,
    terms: Sequence | None = None,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
    max_names: int = MAX_NAMES_GET,
    sources=None,
    min_matched=None,
    shape: str | None = "compact",
) -> dict:
    """The answer as Python objects (used for TSV and by tests; JSON output uses ``render_hits_json``)."""
    shape = parse_shape(shape)
    h = _search(reader, names, terms=terms, include_taxids=include_taxids,
                exclude_taxids=exclude_taxids, alias=alias, max_names=max_names,
                sources=sources, min_matched=min_matched)
    rec = h.store.rec
    if shape == "records":
        results = []
        for o in h.ords:
            aid, taxid, source, date, asm, org, acc, gff = rec[o]
            matched = sorted(h.ann_matched.get(o, []))
            matched_set = set(matched)
            results.append({
                "annotation_id": aid, "organism_name": org, "taxid": taxid,
                "assembly_accession": acc, "assembly_name": asm,
                "source": SOURCE_LABELS.get(source), "release_date": date,
                "n_matched": len(matched), "n_terms": h.n_terms, "matched": matched,
                "alias_only": sorted(h.ann_alias.get(o, [])),
                "missing": [i for i in range(h.n_terms) if i not in matched_set],
                "download_url": download_url(gff),
            })
        return {**h.base, "total": len(results), "facets": h.facets, "results": results}
    rows = []
    for o in h.ords:
        aid, taxid, source, date, asm, org, acc, _gff = rec[o]
        rows.append([aid, taxid, source, sorted(h.ann_matched.get(o, [])), date, asm,
                     sorted(h.ann_alias.get(o, [])), org, acc])
    return {**h.base, "fields": list(HITS_FIELDS), "total": len(rows), "facets": h.facets,
            "rows": rows}


def render_hits_json(
    reader,
    names: str | Sequence[str] | None = None,
    *,
    terms: Sequence | None = None,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
    max_names: int = MAX_NAMES_GET,
    sources=None,
    min_matched=None,
    shape: str | None = "compact",
) -> bytes:
    """Same answer as ``list_name_hits_presence`` as compact JSON bytes, joined from stored pieces."""
    shape = parse_shape(shape)
    h = _search(reader, names, terms=terms, include_taxids=include_taxids,
                exclude_taxids=exclude_taxids, alias=alias, max_names=max_names,
                sources=sources, min_matched=min_matched)
    store, am, aa, n_terms = h.store, h.ann_matched, h.ann_alias, h.n_terms
    total = len(h.ords)
    head: dict = dict(h.base)
    if shape == "records":
        head["total"] = total
        head["facets"] = h.facets
        key = "results"
        parts = []
        for o in h.ords:
            first, last = store.record_pieces(o)
            m = am.get(o, ())
            m_sorted = m if len(m) < 2 else sorted(m)
            if len(m_sorted) == n_terms:
                missing = "[]"
            else:
                have = set(m_sorted)
                missing = _ints([i for i in range(n_terms) if i not in have])
            a = aa.get(o, ())
            parts.append(f'{first}{len(m_sorted)},"n_terms":{n_terms},"matched":{_ints(m_sorted)},'
                         f'"alias_only":{_ints(a if len(a) < 2 else sorted(a))},'
                         f'"missing":{missing}{last}')
    else:
        head["fields"] = list(HITS_FIELDS)
        head["total"] = total
        head["facets"] = h.facets
        key = "rows"
        parts = []
        for o in h.ords:
            first, middle, last = store.compact_pieces(o)
            m = am.get(o, ())
            a = aa.get(o, ())
            parts.append(f"{first}{_ints(m if len(m) < 2 else sorted(m))}{middle}"
                         f"{_ints(a if len(a) < 2 else sorted(a))}{last}")
    prefix = json.dumps(head, ensure_ascii=False, separators=(",", ":"))[:-1]
    return f'{prefix},"{key}":[{",".join(parts)}]}}'.encode("utf-8")


def _count_loci(loci, name_id: int, mode: str, ords: set[int] | None) -> int:
    """Gene rows of one name, optionally restricted to a set of annotation ords (streamed)."""
    sym_only = mode != "include"
    if ords is None:
        sql = "SELECT COUNT(*) FROM locus WHERE name_id=?" + (" AND via=0" if sym_only else "")
        return int(loci.execute(sql, (name_id,)).fetchone()[0])
    n = 0
    for row in loci.execute("SELECT ord, via FROM locus WHERE name_id=?", (name_id,)):
        if row[0] in ords and not (sym_only and row[1] == 1):
            n += 1
    return n


def list_names_by_prefix_presence(
    reader,
    prefix: str,
    *,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
    sources=None,
    loci=None,
) -> dict:
    """Prefix discovery. ``sources`` limits the counted annotations (Ensembl/GenBank/RefSeq);
    ``loci`` (a loci.sqlite connection) adds ``n_loci`` per name: gene rows in the same scope, symbol
    rows only unless alias=include (None when the loci index is not available)."""
    mode = parse_alias_mode(alias)
    source_codes = parse_sources(sources)
    stem = parse_prefix_query(prefix)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")
    conn = reader.conn

    found = conn.execute(
        "SELECT name, name_id, n_sym, n_alias_only FROM name_dict "
        "WHERE name >= ? AND name < ? ORDER BY name LIMIT ?",
        (stem, stem + "￿", MAX_PREFIX_EXPAND + 1),
    ).fetchall()
    if len(found) > MAX_PREFIX_EXPAND:
        raise QueryError(
            f"prefix expands to more than {MAX_PREFIX_EXPAND} accessions",
            code="prefix_too_broad", max=MAX_PREFIX_EXPAND, expanded=len(found), prefix=stem,
        )

    counts: dict[str, tuple[int, int]] = {}
    scope = scope_mask(conn, reader.all_mask, include, exclude)
    if source_codes:
        ph = ",".join("?" for _ in source_codes)
        smask = 0
        for (o,) in conn.execute(f"SELECT ord FROM annotation WHERE source IN ({ph})", source_codes):
            smask |= 1 << int(o)
        scope &= smask
    if include or exclude or source_codes:
        by_id = {int(r["name_id"]): r["name"] for r in found}
        for chunk in _chunks(list(by_id)):
            ph = ",".join("?" for _ in chunk)
            for p in conn.execute(
                f"SELECT name_id, sym, alias FROM posting WHERE name_id IN ({ph})", chunk):
                sym = to_bitset(p["sym"]) & scope
                als = to_bitset(p["alias"]) & scope
                counts[by_id[int(p["name_id"])]] = (sym.bit_count(), (als & ~sym).bit_count())
    else:
        for r in found:
            counts[r["name"]] = (int(r["n_sym"]), int(r["n_alias_only"]))

    ids = {r["name"]: int(r["name_id"]) for r in found}
    scoped_ords = None
    if loci is not None and (include or exclude or source_codes):
        scoped_ords = set(iter_bits(scope))
    results = []
    for name, (n_sym, n_extra) in counts.items():
        n_ann = n_sym + n_extra if mode == "include" else n_sym
        extra = 0 if mode in ("off", "include") else n_extra
        if n_ann == 0:      # alias-only names are listed only with alias=include
            continue
        row = {"name": name, "category": classify_name(name),
               "n_annotations": n_ann, "n_alias_extra": extra}
        if loci is not None:
            row["n_loci"] = _count_loci(loci, ids[name], mode, scoped_ords)
        results.append(row)
    results.sort(key=lambda r: (-r["n_annotations"], -r["n_alias_extra"], r["name"]))
    return {
        "prefix": stem, "alias": mode, "include_taxids": include, "exclude_taxids": exclude,
        "source_filter": [SOURCE_LABELS[c] for c in source_codes],
        "total": len(results), "results": results, "data_version": reader.data_version,
    }
