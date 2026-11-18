"""``/hits`` and ``/names`` over presence.sqlite, with alias as a suggestion (API-DESIGN 3)."""
from __future__ import annotations

from collections.abc import Iterable, Sequence

from helpers.cursor import QueryError
from helpers.names import (
    parse_exact_name,
    parse_name_list,
    parse_prefix_query,
    parse_taxid_list,
)
from helpers.name_category import classify_name
from helpers.pagination import MAX_NAMES_GET, MAX_PREFIX_EXPAND
from helpers.postings import iter_bits, to_bitset
from helpers.presence import scope_mask
from helpers.source_codes import SOURCE_NAME_TO_CODE, SOURCES_LEGEND
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


def list_name_hits_presence(
    reader,
    names: str | Sequence[str],
    *,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
    max_names: int = MAX_NAMES_GET,
    sources=None,
    min_matched=None,
    shape: str | None = "compact",
) -> dict:
    mode = parse_alias_mode(alias)
    shape = parse_shape(shape)
    source_codes = parse_sources(sources)
    raw_names = parse_name_list(names, max_names=max_names)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")

    terms: list[dict] = []
    name_to_term: dict[str, int] = {}
    for raw in raw_names:
        inp, name = parse_exact_name(raw)
        if name in name_to_term:
            continue
        name_to_term[name] = len(terms)
        terms.append({"input": inp, "name": name, "namespaces": [], "n_hits": 0,
                      "n_alias_extra": 0})

    conn = reader.conn
    scope = scope_mask(conn, reader.all_mask, include, exclude)

    ann_matched: dict[int, list[int]] = {}
    ann_alias: dict[int, list[int]] = {}
    for chunk in _chunks(list(name_to_term)):
        ph = ",".join("?" for _ in chunk)
        for row in conn.execute(
            f"SELECT d.name, p.sym, p.alias FROM name_dict d JOIN posting p USING (name_id) "
            f"WHERE d.name IN ({ph})", chunk):
            i = name_to_term[row["name"]]
            sym = to_bitset(row["sym"]) & scope
            als = to_bitset(row["alias"]) & scope
            matched, alias_only, extra = _split_by_mode(sym, als, mode)
            ns = []
            if sym:
                ns.append("symbol")
            if als and mode != "off":
                ns.append("alias")
            terms[i]["namespaces"] = ns
            terms[i]["n_hits"] = matched.bit_count()
            terms[i]["n_alias_extra"] = extra.bit_count()
            for o in iter_bits(matched):
                ann_matched.setdefault(o, []).append(i)
            for o in iter_bits(alias_only):
                ann_alias.setdefault(o, []).append(i)

    n_terms = len(terms)
    needed = parse_min_matched(min_matched, n_terms)
    base = {
        "terms": terms, "sources": dict(HITS_SOURCES), "alias": mode,
        "include_taxids": include, "exclude_taxids": exclude,
        "source_filter": [SOURCE_LABELS[c] for c in source_codes], "min_matched": needed,
        "data_version": reader.data_version,
    }
    ords = sorted(set(ann_matched) | set(ann_alias))
    if needed:
        ords = [o for o in ords if len(ann_matched.get(o, ())) >= needed]
    scored: list[tuple[int, str, dict]] = []
    for chunk in _chunks(ords):
        ph = ",".join("?" for _ in chunk)
        for r in conn.execute(
            "SELECT ord, annotation_id, taxid, organism_name, assembly_accession, "
            f"assembly_name, source, release_date, gff_path FROM annotation WHERE ord IN ({ph})",
            chunk):
            o = int(r["ord"])
            matched = sorted(ann_matched.get(o, []))
            scored.append((len(matched), str(r["annotation_id"]), {
                "annotation_id": str(r["annotation_id"]), "taxid": int(r["taxid"]),
                "source": int(r["source"]), "matched": matched,
                "release_date": _optional_text(r["release_date"]),
                "assembly_name": _optional_text(r["assembly_name"]),
                "alias_only": sorted(ann_alias.get(o, [])),
                "organism_name": _optional_text(r["organism_name"]),
                "assembly_accession": _optional_text(r["assembly_accession"]),
                "gff_path": r["gff_path"],
            }))
    # Counts per source are taken before the sources filter, so a client can show the alternatives.
    facets = {"source": {SOURCE_LABELS[c]: 0 for c in sorted(SOURCE_LABELS)}}
    for _n, _aid, rec in scored:
        label = SOURCE_LABELS.get(rec["source"])
        if label:
            facets["source"][label] += 1
    if source_codes:
        scored = [t for t in scored if t[2]["source"] in source_codes]
    scored.sort(key=lambda t: (-t[0], t[1]))
    recs = [t[2] for t in scored]
    if shape == "records":
        results = []
        for rec in recs:
            matched_set = set(rec["matched"])
            results.append({
                "annotation_id": rec["annotation_id"], "organism_name": rec["organism_name"],
                "taxid": rec["taxid"], "assembly_accession": rec["assembly_accession"],
                "assembly_name": rec["assembly_name"], "source": SOURCE_LABELS.get(rec["source"]),
                "release_date": rec["release_date"], "n_matched": len(rec["matched"]),
                "n_terms": n_terms, "matched": rec["matched"], "alias_only": rec["alias_only"],
                "missing": [i for i in range(n_terms) if i not in matched_set],
                "download_url": download_url(rec["gff_path"]),
            })
        return {**base, "total": len(results), "facets": facets, "results": results}
    rows = [[rec["annotation_id"], rec["taxid"], rec["source"], rec["matched"],
             rec["release_date"], rec["assembly_name"], rec["alias_only"], rec["organism_name"],
             rec["assembly_accession"]] for rec in recs]
    return {**base, "fields": list(HITS_FIELDS), "total": len(rows), "facets": facets, "rows": rows}


def list_names_by_prefix_presence(
    reader,
    prefix: str,
    *,
    include_taxids=None,
    exclude_taxids=None,
    alias: str | None = "suggest",
) -> dict:
    mode = parse_alias_mode(alias)
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
    if include or exclude:
        scope = scope_mask(conn, reader.all_mask, include, exclude)
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

    results = []
    for name, (n_sym, n_extra) in counts.items():
        n_ann = n_sym + n_extra if mode == "include" else n_sym
        extra = 0 if mode in ("off", "include") else n_extra
        if n_ann == 0:      # alias-only names are listed only with alias=include
            continue
        results.append({"name": name, "category": classify_name(name),
                        "n_annotations": n_ann, "n_alias_extra": extra})
    results.sort(key=lambda r: (-r["n_annotations"], -r["n_alias_extra"], r["name"]))
    return {
        "prefix": stem, "alias": mode, "include_taxids": include, "exclude_taxids": exclude,
        "total": len(results), "results": results, "data_version": reader.data_version,
    }
