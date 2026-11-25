"""``/distributions``: how long a gene is, per annotation, grouped by clade.

Two views of the same numbers, in the same compact-table shape as ``/hits``:

* ``summarise``: one row per (name, column), a column being a clade the client chose;
* ``points``: one row per annotation for one name in one clade.

Length is the gene span ``end - start + 1`` as stored in ``loci.sqlite``. An annotation counts once
per name: when several genes match (copies), the value is their mean length and ``copies`` says how
many. Memory is flat: the rows of one name are streamed, folded into per-column accumulators and
dropped before the next name; the largest thing held is one name's per-annotation values
(at most one float per annotation).
"""
from __future__ import annotations

import math
import sqlite3
from collections.abc import Iterator, Sequence
from itertools import groupby
from operator import itemgetter

from helpers.cursor import QueryError
from helpers.names import parse_taxid_list
from helpers.pagination import MAX_TAXIDS_FILTER
from helpers.postings import iter_bits
from helpers.presence import scope_mask
from services.presence_loci import LociQuery

MAX_COLUMNS = 50

# Fixed log-size buckets, shared by every cell so the client can add them up or compare them.
BIN_COUNT = 40
BIN_MIN_BP = 100
BIN_MAX_BP = 10_000_000
_LOG_MIN = math.log10(BIN_MIN_BP)
_LOG_SPAN = math.log10(BIN_MAX_BP) - _LOG_MIN

BINS = {"scale": "log10", "min_bp": BIN_MIN_BP, "max_bp": BIN_MAX_BP, "n": BIN_COUNT}

# ``p5`` and ``p95`` come after ``hist`` so clients written for the first version keep working;
# they are the whiskers of a box plot that outliers do not stretch.
SUMMARY_FIELDS = ("name", "column", "n_with_gene", "n_multi_copy", "mean", "min", "p25", "median",
                  "p75", "max", "hist", "p5", "p95")
POINT_FIELDS = ("annotation_id", "taxid", "length", "copies", "via")


def bin_of(length: float) -> int:
    """Bucket of a length: lengths below ``BIN_MIN_BP`` or above ``BIN_MAX_BP`` go to the end buckets."""
    if length <= BIN_MIN_BP:
        return 0
    if length >= BIN_MAX_BP:
        return BIN_COUNT - 1
    return min(BIN_COUNT - 1, int((math.log10(length) - _LOG_MIN) / _LOG_SPAN * BIN_COUNT))


def quantile(sorted_values: Sequence[float], p: float) -> float:
    """``p`` in [0, 1], linear interpolation between neighbours (the usual "type 7")."""
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    pos = p * (n - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 >= n:
        return float(sorted_values[-1])
    return sorted_values[lo] + (sorted_values[lo + 1] - sorted_values[lo]) * frac


def parse_columns(raw) -> list[dict]:
    """``[{"id": "mammals", "taxids": [40674]}, ...]`` -> ``[{"id", "taxids"}]``.

    A column without taxids covers every annotation in scope. No column at all means one column
    ``all``.
    """
    if raw is None:
        return [{"id": "all", "taxids": []}]
    if not isinstance(raw, list) or not raw:
        raise QueryError("columns must be a non-empty list", code="invalid_columns")
    if len(raw) > MAX_COLUMNS:
        raise QueryError(f"columns accepts at most {MAX_COLUMNS} values", code="too_many_columns",
                         max=MAX_COLUMNS, received=len(raw))
    out: list[dict] = []
    seen: set[str] = set()
    for i, item in enumerate(raw):
        if hasattr(item, "model_dump"):
            item = item.model_dump()
        if not isinstance(item, dict):
            raise QueryError("each column must be an object with taxids", code="invalid_columns",
                             index=i)
        taxids = parse_taxid_list(item.get("taxids"), field="taxids", max_taxids=MAX_TAXIDS_FILTER)
        cid = item.get("id")
        cid = str(cid).strip() if cid not in (None, "") else (",".join(map(str, taxids)) or "all")
        if cid in seen:
            raise QueryError(f"column id {cid!r} is used twice", code="duplicate_column", id=cid)
        seen.add(cid)
        out.append({"id": cid, "taxids": taxids})
    return out


def names_of(q: LociQuery) -> list[dict]:
    """The names a request expands to, in term order: ``[{"name", "name_id", "term"}]``."""
    out, seen = [], set()
    for term in q.terms:
        for name, name_id in zip(term["names"], term["name_ids"]):
            if name_id not in seen:
                seen.add(name_id)
                out.append({"name": name, "name_id": name_id, "term": term["index"]})
    return out


def annotation_values(loci: sqlite3.Connection, name_id: int, scope: set[int],
                      alias_mode: str) -> Iterator[tuple[int, float, int, int]]:
    """``(ord, mean_length, copies, via)`` per annotation of ``scope`` that has this name.

    ``via`` is 0 when the value comes from genes whose symbol is the name and 1 when it comes
    from genes that only carry it as an alias. Symbol genes win: alias genes are used only when the
    annotation has no symbol gene, and never with ``alias=off``.
    """
    rows = loci.execute("SELECT ord, start, end, via FROM locus WHERE name_id=? "
                        "ORDER BY ord, local_id", (name_id,))
    for ord_, grp in groupby(rows, key=itemgetter(0)):
        if ord_ not in scope:
            continue
        symbol, alias = [], []
        for _o, start, end, via in grp:
            (alias if via == 1 else symbol).append(end - start + 1)
        if symbol:
            yield ord_, sum(symbol) / len(symbol), len(symbol), 0
        elif alias and alias_mode != "off":
            yield ord_, sum(alias) / len(alias), len(alias), 1


def column_masks(pres: sqlite3.Connection, all_mask: int, columns: list[dict],
                 scope: set[int]) -> list[set[int]]:
    """Annotations (ordinals) of each column, limited to the request's scope."""
    out = []
    for col in columns:
        mask = scope_mask(pres, all_mask, col["taxids"], [])
        out.append({o for o in iter_bits(mask) if o in scope})
    return out


def summarise(loci: sqlite3.Connection, q: LociQuery, names: list[dict],
              col_sets: list[set[int]]) -> list[list]:
    """Rows ``[name_i, column_i, n_with_gene, n_multi_copy, mean, min, p25, median, p75, max, hist, p5, p95]``.

    Cells without the gene are left out. Columns may overlap; an annotation counts in each column
    that holds it.
    """
    cols_of: dict[int, list[int]] = {}
    for ci, members in enumerate(col_sets):
        for o in members:
            cols_of.setdefault(o, []).append(ci)
    rows: list[list] = []
    for ni, item in enumerate(names):
        values: list[list[float]] = [[] for _ in col_sets]
        multi = [0] * len(col_sets)
        hist = [[0] * BIN_COUNT for _ in col_sets]
        for ord_, mean_len, copies, _via in annotation_values(loci, item["name_id"], q.scope, q.mode):
            for ci in cols_of.get(ord_, ()):
                values[ci].append(mean_len)
                hist[ci][bin_of(mean_len)] += 1
                if copies > 1:
                    multi[ci] += 1
        for ci, vals in enumerate(values):
            if not vals:
                continue
            vals.sort()
            rows.append([ni, ci, len(vals), multi[ci], round(sum(vals) / len(vals), 1),
                         round(vals[0], 1), round(quantile(vals, 0.25), 1),
                         round(quantile(vals, 0.5), 1), round(quantile(vals, 0.75), 1),
                         round(vals[-1], 1), hist[ci], round(quantile(vals, 0.05), 1),
                         round(quantile(vals, 0.95), 1)])
    return rows


def points(loci: sqlite3.Connection, q: LociQuery, name_id: int,
           ann_info: dict[int, tuple], ann_paths: dict[int, tuple[int, str]]) -> list[list]:
    """Rows ``[annotation_id, taxid, length, copies, via]`` of one name, in annotation order."""
    rows = []
    for ord_, mean_len, copies, via in annotation_values(loci, name_id, q.scope, q.mode):
        rows.append([ann_info[ord_][0], ann_paths[ord_][0], round(mean_len, 1), copies, via])
    return rows
