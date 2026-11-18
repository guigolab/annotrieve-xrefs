"""``/loci``: where the genes are, from loci.sqlite (API-DESIGN 4.3).

Rows come out in a deterministic order (annotation ordinal, then seqid, start, end, local_id),
one row per gene even when several terms hit it. Memory is bounded by one annotation's rows.
"""
from __future__ import annotations

import base64
import binascii
import hashlib
import heapq
import json
import sqlite3
import threading
from collections import OrderedDict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from itertools import groupby
from operator import itemgetter
from urllib.parse import quote

from helpers.cursor import QueryError
from helpers.names import parse_exact_name, parse_name_list, parse_taxid_list
from helpers.pagination import MAX_NAMES_GET
from helpers.postings import iter_bits, to_bitset
from helpers.presence import scope_mask
from helpers.source_codes import SOURCES_LEGEND
from services.terms import parse_terms, resolve_terms
from services.presence_hits import _chunks, parse_alias_mode
from settings import settings

MAX_FLANK_BP = 1_000_000
STRAND = {1: "+", -1: "-", 0: "."}


@dataclass
class LociQuery:
    mode: str
    case_sensitive: bool
    flank_bp: int
    include: list[int]
    exclude: list[int]
    annotation_ids: list[str]
    terms: list[dict]
    name_ids: list[list[int]]          # per term: the index names it expands to ([] = not found)
    scope: set[int]
    data_version: str
    scoped: bool = False
    key: str = field(default="")
    id_name: dict[int, str] = field(default_factory=dict)


def parse_flank(value) -> int:
    try:
        flank = int(value or 0)
    except (TypeError, ValueError):
        raise QueryError("flank_bp must be an integer", code="invalid_flank", value=value) from None
    if flank < 0 or flank > MAX_FLANK_BP:
        raise QueryError(f"flank_bp must be between 0 and {MAX_FLANK_BP}",
                         code="invalid_flank", value=flank)
    return flank


def parse_annotation_ids(raw) -> list[str]:
    if raw is None:
        return []
    pieces = [raw] if isinstance(raw, str) else list(raw)
    out: list[str] = []
    for piece in pieces:
        for part in str(piece).split(","):
            part = part.strip()
            if part and part not in out:
                out.append(part)
    if len(out) > 5000:
        raise QueryError("annotation_ids accepts at most 5000 values",
                         code="too_many_annotation_ids", max=5000, received=len(out))
    return out


def build_query(pres: sqlite3.Connection, all_mask: int, data_version: str, names=None,
                *, terms=None, include_taxids=None, exclude_taxids=None, alias="suggest",
                annotation_ids=None, case_sensitive=False, flank_bp=0,
                max_names: int = MAX_NAMES_GET) -> LociQuery:
    mode = parse_alias_mode(alias)
    raws = parse_terms(names, terms, max_terms=max_names)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")
    ann_ids = parse_annotation_ids(annotation_ids)
    flank = parse_flank(flank_bp)

    resolved = resolve_terms(pres, raws, mode)
    name_ids = [t["name_ids"] for t in resolved]
    id_name = {nid: n for t in resolved for n, nid in zip(t["names"], t["name_ids"])}

    mask = scope_mask(pres, all_mask, include, exclude)
    if ann_ids:
        allowed = 0
        for chunk in _chunks(ann_ids):
            ph = ",".join("?" for _ in chunk)
            for (o,) in pres.execute(f"SELECT ord FROM annotation WHERE annotation_id IN ({ph})", chunk):
                allowed |= 1 << int(o)
        mask &= allowed
    key = hashlib.sha1(json.dumps(
        [[[t["kind"], t["name"], t["names"]] for t in resolved], include, exclude, ann_ids, mode,
         case_sensitive, [t["typed"] for t in resolved] if case_sensitive else None, data_version],
        sort_keys=True).encode()).hexdigest()[:12]
    return LociQuery(mode, bool(case_sensitive), flank, include, exclude, ann_ids, resolved,
                     name_ids, set(iter_bits(mask)), data_version,
                     scoped=bool(include or exclude or ann_ids), key=key, id_name=id_name)


def encode_cursor(q: LociQuery, ord_: int, n: int) -> str:
    raw = json.dumps({"o": ord_, "n": n, "q": q.key}).encode()
    return base64.urlsafe_b64encode(raw).decode().rstrip("=")


def decode_cursor(q: LociQuery, cursor: str | None) -> tuple[int, int] | None:
    if not cursor:
        return None
    try:
        data = json.loads(base64.urlsafe_b64decode(cursor + "=" * (-len(cursor) % 4)))
        ord_, n, key = int(data["o"]), int(data["n"]), data["q"]
    except (ValueError, KeyError, TypeError, binascii.Error):
        raise QueryError("invalid cursor", code="invalid_cursor") from None
    if key != q.key:
        raise QueryError("cursor does not match this query", code="cursor_filter_mismatch")
    return ord_, n


def estimate(pres: sqlite3.Connection, loci: sqlite3.Connection, q: LociQuery) -> dict:
    """``n_annotations`` (exact, from the postings) and ``n_loci_upper_bound`` (unscoped)."""
    ids = sorted({i for term in q.name_ids for i in term})
    mask = 0
    for chunk in _chunks(ids):
        ph = ",".join("?" for _ in chunk)
        for r in pres.execute(f"SELECT sym, alias FROM posting WHERE name_id IN ({ph})", chunk):
            m = to_bitset(r["sym"])
            if q.mode != "off":
                m |= to_bitset(r["alias"])
            mask |= m
    scope = 0
    for o in q.scope:
        scope |= 1 << o
    bound = 0
    for chunk in _chunks(ids):
        ph = ",".join("?" for _ in chunk)
        bound += loci.execute(
            f"SELECT COALESCE(SUM(n_loci), 0) FROM name_loci WHERE name_id IN ({ph})", chunk
        ).fetchone()[0]
    return {"n_annotations": (mask & scope).bit_count(), "n_loci_upper_bound": int(bound)}


def gff_url(annotation_id: str, seqid: str, start: int, end: int, flank: int) -> str:
    """Annotrieve slice link. pysam is 0-based half-open: GFF start S..E -> start=S-1, end=E."""
    base = settings.ANNOTRIEVE_API_BASE.rstrip("/")
    return (f"{base}/annotations/{quote(annotation_id, safe='')}/gff?region={quote(seqid, safe='')}"
            f"&start={max(0, start - 1 - flank)}&end={end + flank}")


def _term_rows(loci: sqlite3.Connection, term_i: int, name_ids: list[int], min_ord: int,
               q: LociQuery) -> Iterator[tuple]:
    """Rows of one term, ordered by (ord, local_id). A gene hit by several of the term's names
    appears once per name; the caller merges them."""
    ph = ",".join("?" for _ in name_ids)
    cur = loci.execute(
        "SELECT ord, local_id, seqid_id, start, end, strand, feature_type_id, biotype_id, via, "
        f"gene_id, primary_name, name_id FROM locus WHERE name_id IN ({ph}) AND ord>=? "
        "ORDER BY ord, local_id", (*name_ids, min_ord))
    term = q.terms[term_i]
    typed = term["typed"]
    whole = term["kind"] == "exact" and not term.get("from_prefix")
    for r in cur:
        if r["ord"] not in q.scope:
            continue
        via = r["via"]
        if via == 1 and q.mode == "off":
            continue
        if q.case_sensitive:
            pn = r["primary_name"] or ""
            if via == 1 or not (pn == typed if whole else pn.startswith(typed)):
                continue
        yield (r["ord"], r["local_id"], term_i, r)


_SEQID_CACHE: "OrderedDict[tuple[str, int, int], str]" = OrderedDict()
_SEQID_LOCK = threading.Lock()
SEQID_CACHE_ENTRIES = 50_000       # one short string each (about 10 MB at most)


def _seqid(loci: sqlite3.Connection, data_version: str, ord_: int, seqid_id: int) -> str:
    """One contig name by key. Whole per-annotation tables are not loaded: a scaffold-level
    assembly can have tens of thousands of contigs and a query touches only a few of them."""
    key = (data_version, ord_, seqid_id)
    with _SEQID_LOCK:
        hit = _SEQID_CACHE.get(key)
        if hit is not None:
            _SEQID_CACHE.move_to_end(key)
            return hit
    row = loci.execute("SELECT seqid FROM seqid WHERE ord=? AND seqid_id=?",
                       (ord_, seqid_id)).fetchone()
    value = row[0] if row else ""
    with _SEQID_LOCK:
        _SEQID_CACHE[key] = value
        while len(_SEQID_CACHE) > SEQID_CACHE_ENTRIES:
            _SEQID_CACHE.popitem(last=False)
    return value


def iter_loci(pres: sqlite3.Connection, loci: sqlite3.Connection, q: LociQuery,
              after: tuple[int, int] | None = None,
              ann_info: dict[int, tuple] | None = None) -> Iterator[dict]:
    """Yield locus dicts in deterministic order; ``after`` = (ord, rows already sent in it).

    ``ann_info`` (ord -> annotation_id, organism_name, assembly_accession, source) avoids one
    query per annotation; without it the presence connection is queried.
    """
    min_ord, skip = (after or (0, 0))
    ft = dict(loci.execute("SELECT id, name FROM dict_feature_type"))
    bt = dict(loci.execute("SELECT id, name FROM dict_biotype"))
    streams = [_term_rows(loci, i, ids, min_ord, q)
               for i, ids in enumerate(q.name_ids) if ids]
    merged = heapq.merge(*streams, key=itemgetter(0, 1))
    for ord_, grp in groupby(merged, key=itemgetter(0)):
        genes: dict[int, dict] = {}
        for _o, local_id, term_i, r in grp:
            g = genes.get(local_id)
            if g is None:
                g = genes[local_id] = {"row": r, "terms": [], "names": set(), "via": 1}
            g["terms"].append(term_i)
            g["names"].add(q.id_name[r["name_id"]])
            if r["via"] == 0:
                g["via"] = 0
        if ann_info is not None:
            info = ann_info[ord_]
        else:
            r0 = pres.execute("SELECT annotation_id, organism_name, assembly_accession, source "
                              "FROM annotation WHERE ord=?", (ord_,)).fetchone()
            info = (r0["annotation_id"], r0["organism_name"], r0["assembly_accession"],
                    int(r0["source"]))
        ann_id, organism, accession, source = info
        out = []
        for local_id, g in genes.items():
            r = g["row"]
            terms = sorted(set(g["terms"]))
            matched = sorted(g["names"])
            pname = r["primary_name"]
            symbol = pname if pname and pname.casefold() in g["names"] else matched[0]
            seqid = _seqid(loci, q.data_version, ord_, r["seqid_id"])
            out.append({
                "annotation_id": ann_id, "ord": ord_,
                "organism_name": organism,
                "assembly_accession": accession,
                "source": SOURCES_LEGEND.get(str(source), ""),
                "locus_ref": f"{ann_id}:{local_id}", "local_id": local_id,
                "gene_id": r["gene_id"], "symbol": symbol, "primary_name": pname,
                "terms": terms, "matched_names": matched,
                "matched_via": "symbol" if g["via"] == 0 else "alias",
                "seqid": seqid, "start": r["start"], "end": r["end"],
                "strand": STRAND.get(r["strand"], "."),
                "feature_type": ft.get(r["feature_type_id"]),
                "biotype": bt.get(r["biotype_id"]),
                "gff_url": gff_url(ann_id, seqid, r["start"], r["end"], q.flank_bp),
            })
        out.sort(key=lambda d: (d["seqid"], d["start"], d["end"], d["local_id"]))
        start_at = skip if ord_ == min_ord else 0
        for row in out[start_at:]:
            yield row
