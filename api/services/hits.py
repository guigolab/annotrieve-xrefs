"""Gene-name presence search (``/hits``)."""
from __future__ import annotations

import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Sequence

from helpers.cursor import QueryError
from helpers.match import MATCH_ALL, MATCH_ANY, MatchMode, parse_match
from helpers.names import (
    MAX_PREFIX_EXPAND,
    parse_name_list,
    parse_taxid_list,
    split_name_term,
)
from helpers.pagination import MAX_NAMES_GET
from helpers.source_codes import SOURCES_LEGEND, source_code

HITS_FIELDS = [
    "annotation_id",
    "taxid",
    "source",
    "n_matched",
    "matched",
    "release_date",
    "assembly_name",
]

HITS_SOURCES: dict[str, dict[str, str]] = {
    code: {"label": label} for code, label in SOURCES_LEGEND.items()
}

# Stay under common SQLite SQLITE_MAX_VARIABLE_NUMBER (999 on some builds).
_SQLITE_MAX_VARS = 900


def _chunks(items: list, size: int) -> Iterable[list]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1",
        (name,),
    ).fetchone()
    return row is not None


def _expand_prefix(conn: sqlite3.Connection, prefix: str) -> list[str]:
    """
    Distinct symbol∪alias accessions with accession prefix, ranked by coverage.

    Raises ``prefix_too_broad`` when more than ``MAX_PREFIX_EXPAND`` distinct
    accessions match.
    """
    rows = conn.execute(
        """
        SELECT accession, MAX(n_annotations) AS n_ann
        FROM xref_meta
        WHERE namespace IN ('symbol', 'alias')
          AND accession >= ?
          AND accession < ?
        GROUP BY accession
        ORDER BY n_ann DESC, accession ASC
        LIMIT ?
        """,
        (prefix, prefix + "\uffff", MAX_PREFIX_EXPAND + 1),
    ).fetchall()
    if len(rows) > MAX_PREFIX_EXPAND:
        raise QueryError(
            f"prefix expands to more than {MAX_PREFIX_EXPAND} accessions",
            code="prefix_too_broad",
            max=MAX_PREFIX_EXPAND,
            expanded=len(rows),
            prefix=prefix,
        )
    return [str(r["accession"]) for r in rows]


def _namespaces_for_accessions(
    conn: sqlite3.Connection, accessions: list[str]
) -> dict[str, list[str]]:
    """
    Map accession → sorted namespaces present in xref_file (symbol/alias only).
    """
    out: dict[str, list[str]] = {acc: [] for acc in accessions}
    if not accessions:
        return out
    # One bind per accession; chunk under limit.
    for chunk in _chunks(accessions, _SQLITE_MAX_VARS):
        placeholders = ",".join("?" for _ in chunk)
        rows = conn.execute(
            f"""
            SELECT DISTINCT namespace, accession
            FROM xref_file
            WHERE accession IN ({placeholders})
              AND namespace IN ('symbol', 'alias')
            ORDER BY accession, namespace
            """,
            chunk,
        ).fetchall()
        for row in rows:
            acc = str(row["accession"])
            if acc in out:
                out[acc].append(str(row["namespace"]))
    return out


def _postings_for_keys(
    conn: sqlite3.Connection,
    pairs: list[tuple[str, str]],
) -> dict[tuple[str, str], set[int]]:
    """Map (namespace, accession) → annotation_key set (chunked OR queries)."""
    out: dict[tuple[str, str], set[int]] = {k: set() for k in pairs}
    if not pairs:
        return out
    # Each pair uses 2 binds.
    pair_chunk = max(1, _SQLITE_MAX_VARS // 2)
    for chunk in _chunks(pairs, pair_chunk):
        preds = " OR ".join(
            "(namespace = ? AND accession = ?)" for _ in chunk
        )
        params: list = []
        for ns, acc in chunk:
            params.extend([ns, acc])
        rows = conn.execute(
            f"""
            SELECT namespace, accession, annotation_key
            FROM xref_file
            WHERE {preds}
            """,
            params,
        ).fetchall()
        for row in rows:
            key = (str(row["namespace"]), str(row["accession"]))
            if key in out:
                out[key].add(int(row["annotation_key"]))
    return out


def _filter_keys_by_tax(
    conn: sqlite3.Connection,
    annotation_keys: set[int],
    include_taxids: list[int],
    exclude_taxids: list[int],
) -> set[int]:
    if not annotation_keys:
        return set()
    if not include_taxids and not exclude_taxids:
        return set(annotation_keys)

    key_list = list(annotation_keys)

    if include_taxids:
        allowed: set[int] = set()
        inc_ph = ",".join("?" for _ in include_taxids)
        key_chunk = max(1, _SQLITE_MAX_VARS - len(include_taxids))
        for chunk in _chunks(key_list, key_chunk):
            placeholders = ",".join("?" for _ in chunk)
            rows = conn.execute(
                f"""
                SELECT DISTINCT annotation_key
                FROM annotation_lineage
                WHERE annotation_key IN ({placeholders})
                  AND taxid IN ({inc_ph})
                """,
                [*chunk, *include_taxids],
            ).fetchall()
            allowed.update(int(r[0]) for r in rows)
    else:
        allowed = set(annotation_keys)

    if exclude_taxids and allowed:
        banned: set[int] = set()
        exc_ph = ",".join("?" for _ in exclude_taxids)
        key_chunk = max(1, _SQLITE_MAX_VARS - len(exclude_taxids))
        allowed_list = list(allowed)
        for chunk in _chunks(allowed_list, key_chunk):
            placeholders = ",".join("?" for _ in chunk)
            rows = conn.execute(
                f"""
                SELECT DISTINCT annotation_key
                FROM annotation_lineage
                WHERE annotation_key IN ({placeholders})
                  AND taxid IN ({exc_ph})
                """,
                [*chunk, *exclude_taxids],
            ).fetchall()
            banned.update(int(r[0]) for r in rows)
        allowed -= banned

    return allowed


def _fetch_annotations(
    conn: sqlite3.Connection, annotation_keys: list[int]
) -> dict[int, sqlite3.Row]:
    by_key: dict[int, sqlite3.Row] = {}
    if not annotation_keys:
        return by_key
    for chunk in _chunks(annotation_keys, _SQLITE_MAX_VARS):
        placeholders = ",".join("?" for _ in chunk)
        rows = conn.execute(
            f"""
            SELECT annotation_key, annotation_id, taxid, source_database,
                   release_date, assembly_name
            FROM annotation
            WHERE annotation_key IN ({placeholders})
            """,
            chunk,
        ).fetchall()
        for row in rows:
            by_key[int(row["annotation_key"])] = row
    return by_key


def _optional_text(raw) -> str | None:
    if raw is None:
        return None
    text = str(raw).strip()
    return text or None


def _envelope(
    *,
    terms: list[dict],
    keys: list[dict],
    include: list[int],
    exclude: list[int],
    rows: list,
) -> dict:
    return {
        "terms": terms,
        "keys": keys,
        "sources": dict(HITS_SOURCES),
        "include_taxids": include,
        "exclude_taxids": exclude,
        "total": len(rows),
        "fields": list(HITS_FIELDS),
        "rows": rows,
    }


def list_name_hits(
    conn: sqlite3.Connection,
    names: str | Sequence[str],
    *,
    include_taxids: str | Sequence[str] | Sequence[int] | None = None,
    exclude_taxids: str | Sequence[str] | Sequence[int] | None = None,
    match: str | None = MATCH_ANY,
    max_names: int = MAX_NAMES_GET,
) -> dict:
    """
    Expand gene names against symbol∪alias and return presence hits.

    See plan envelope: terms, keys, sources, taxid echoes, fields, rows.
    """
    match_mode: MatchMode = parse_match(match)
    raw_names = parse_name_list(names, max_names=max_names)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")

    if not _table_exists(conn, "xref_file") or not _table_exists(
        conn, "xref_meta"
    ):
        terms_out = []
        for raw in raw_names:
            _inp, name_form, is_prefix = split_name_term(raw)
            term: dict = {
                "input": _inp,
                "name": name_form,
                "matched_namespaces": [],
                "n_hits": 0,
            }
            if is_prefix:
                term["expanded"] = []
            terms_out.append(term)
        return _envelope(
            terms=terms_out, keys=[], include=include, exclude=exclude, rows=[]
        )

    # --- expand terms → candidate accessions ---
    term_specs: list[dict] = []
    all_candidates: list[str] = []
    for raw in raw_names:
        inp, name_form, is_prefix = split_name_term(raw)
        if is_prefix:
            stem = name_form[:-1]
            expanded = _expand_prefix(conn, stem)
            candidates = expanded
        else:
            expanded = []
            candidates = [name_form]
        term_specs.append(
            {
                "input": inp,
                "name": name_form,
                "is_prefix": is_prefix,
                "expanded": expanded,
                "candidates": candidates,
            }
        )
        all_candidates.extend(candidates)

    # Batch accession → namespaces (dedupe query input, preserve per-term loops)
    unique_accs = list(dict.fromkeys(all_candidates))
    acc_to_ns = _namespaces_for_accessions(conn, unique_accs)

    # --- build keys (term_idx, name, ns) present in xref_file ---
    keys: list[dict] = []
    key_ns_acc: list[tuple[str, str]] = []
    term_key_indices: list[list[int]] = [[] for _ in term_specs]

    for term_i, spec in enumerate(term_specs):
        for acc in spec["candidates"]:
            for ns in acc_to_ns.get(acc, []):
                idx = len(keys)
                keys.append({"term": term_i, "name": acc, "ns": ns})
                key_ns_acc.append((ns, acc))
                term_key_indices[term_i].append(idx)

    unique_pairs = list(dict.fromkeys(key_ns_acc))
    postings = _postings_for_keys(conn, unique_pairs)

    all_ann_keys: set[int] = set()
    for s in postings.values():
        all_ann_keys |= s

    allowed = _filter_keys_by_tax(conn, all_ann_keys, include, exclude)

    # Intersect in place; drop unfiltered map reference for GC
    for pair, ann_set in postings.items():
        ann_set &= allowed
    filtered_postings = postings
    del all_ann_keys
    del allowed

    key_ann_sets: list[set[int]] = [
        filtered_postings.get((ns, acc), set()) for ns, acc in key_ns_acc
    ]

    # terms[] — matched_namespaces only from keys with post-tax hits
    terms_out: list[dict] = []
    terms_with_keys: list[int] = []
    for term_i, spec in enumerate(term_specs):
        ns_set: set[str] = set()
        ann_union: set[int] = set()
        for ki in term_key_indices[term_i]:
            if key_ann_sets[ki]:
                ns_set.add(keys[ki]["ns"])
                ann_union |= key_ann_sets[ki]
        term = {
            "input": spec["input"],
            "name": spec["name"],
            "matched_namespaces": sorted(ns_set),
            "n_hits": len(ann_union),
        }
        if spec["is_prefix"]:
            term["expanded"] = list(spec["expanded"])
        terms_out.append(term)
        if term_key_indices[term_i]:
            terms_with_keys.append(term_i)

    if match_mode == MATCH_ALL:
        for term_i in range(len(term_specs)):
            if not term_key_indices[term_i]:
                return _envelope(
                    terms=terms_out,
                    keys=keys,
                    include=include,
                    exclude=exclude,
                    rows=[],
                )

    ann_to_key_idxs: dict[int, list[int]] = defaultdict(list)
    for ki, ann_set in enumerate(key_ann_sets):
        for ak in ann_set:
            ann_to_key_idxs[ak].append(ki)

    row_ann_keys: list[int] = []
    for ak, key_idxs in ann_to_key_idxs.items():
        terms_hit = {keys[ki]["term"] for ki in key_idxs}
        if match_mode == MATCH_ALL:
            if not all(t in terms_hit for t in terms_with_keys):
                continue
        row_ann_keys.append(ak)

    if not row_ann_keys:
        return _envelope(
            terms=terms_out,
            keys=keys,
            include=include,
            exclude=exclude,
            rows=[],
        )

    by_key = _fetch_annotations(conn, row_ann_keys)

    scored: list[tuple[int, str, list]] = []
    for ak in row_ann_keys:
        row = by_key.get(ak)
        if row is None:
            continue
        key_idxs = sorted(ann_to_key_idxs[ak])
        n_matched = len({keys[ki]["term"] for ki in key_idxs})
        scored.append(
            (
                n_matched,
                str(row["annotation_id"]),
                [
                    str(row["annotation_id"]),
                    int(row["taxid"]),
                    source_code(row["source_database"]),
                    n_matched,
                    key_idxs,
                    _optional_text(row["release_date"]),
                    _optional_text(row["assembly_name"]),
                ],
            )
        )

    scored.sort(key=lambda t: (-t[0], t[1]))
    rows = [t[2] for t in scored]

    return _envelope(
        terms=terms_out,
        keys=keys,
        include=include,
        exclude=exclude,
        rows=rows,
    )
