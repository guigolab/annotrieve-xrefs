"""Gene-name presence search (``/hits``) — exact names only."""
from __future__ import annotations

import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Sequence

from helpers.cursor import QueryError
from helpers.names import parse_exact_name, parse_name_list, parse_taxid_list
from helpers.pagination import MAX_NAMES_GET
from helpers.source_codes import SOURCES_LEGEND, source_code

HITS_FIELDS = [
    "annotation_id",
    "taxid",
    "source",
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


def _postings_for_accessions(
    conn: sqlite3.Connection,
    accessions: list[str],
) -> list[tuple[str, str, int]]:
    """Return ``(namespace, accession, annotation_key)`` for symbol/alias."""
    out: list[tuple[str, str, int]] = []
    if not accessions:
        return out
    for chunk in _chunks(accessions, _SQLITE_MAX_VARS):
        placeholders = ",".join("?" for _ in chunk)
        rows = conn.execute(
            f"""
            SELECT namespace, accession, annotation_key
            FROM xref_file
            WHERE accession IN ({placeholders})
              AND namespace IN ('symbol', 'alias')
            """,
            chunk,
        ).fetchall()
        for row in rows:
            out.append(
                (
                    str(row["namespace"]),
                    str(row["accession"]),
                    int(row["annotation_key"]),
                )
            )
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
    include: list[int],
    exclude: list[int],
    rows: list,
) -> dict:
    return {
        "terms": terms,
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
    max_names: int = MAX_NAMES_GET,
) -> dict:
    """
    Exact gene-name presence search against symbol∪alias.

    Envelope: ``terms`` (deduped), ``sources``, taxid echoes, ``fields``,
    ``rows`` whose ``matched`` lists term indices.
    """
    raw_names = parse_name_list(names, max_names=max_names)
    include = parse_taxid_list(include_taxids, field="include_taxids")
    exclude = parse_taxid_list(exclude_taxids, field="exclude_taxids")

    # Dedupe by casefolded name; keep first input spelling.
    terms_out: list[dict] = []
    name_to_term: dict[str, int] = {}
    for raw in raw_names:
        inp, name = parse_exact_name(raw)
        if name in name_to_term:
            continue
        name_to_term[name] = len(terms_out)
        terms_out.append(
            {
                "input": inp,
                "name": name,
                "namespaces": [],
                "n_hits": 0,
            }
        )

    if not _table_exists(conn, "xref_file"):
        return _envelope(
            terms=terms_out, include=include, exclude=exclude, rows=[]
        )

    accessions = list(name_to_term.keys())
    postings = _postings_for_accessions(conn, accessions)

    all_ann_keys = {ak for _ns, _acc, ak in postings}
    allowed = _filter_keys_by_tax(conn, all_ann_keys, include, exclude)

    candidate_keys = {ak for _ns, _acc, ak in postings if ak in allowed}
    if not candidate_keys:
        return _envelope(
            terms=terms_out, include=include, exclude=exclude, rows=[]
        )

    by_key = _fetch_annotations(conn, list(candidate_keys))

    # term_i → (namespaces with hits, annotation_key set); only real annotations.
    term_ns: list[set[str]] = [set() for _ in terms_out]
    term_anns: list[set[int]] = [set() for _ in terms_out]
    ann_to_terms: dict[int, set[int]] = defaultdict(set)

    for ns, acc, ak in postings:
        if ak not in by_key:
            continue
        term_i = name_to_term.get(acc)
        if term_i is None:
            continue
        term_ns[term_i].add(ns)
        term_anns[term_i].add(ak)
        ann_to_terms[ak].add(term_i)

    for term_i, term in enumerate(terms_out):
        term["namespaces"] = sorted(term_ns[term_i])
        term["n_hits"] = len(term_anns[term_i])

    if not ann_to_terms:
        return _envelope(
            terms=terms_out, include=include, exclude=exclude, rows=[]
        )

    scored: list[tuple[int, str, list]] = []
    for ak, matched_set in ann_to_terms.items():
        row = by_key[ak]
        matched = sorted(matched_set)
        n_matched = len(matched)
        scored.append(
            (
                n_matched,
                str(row["annotation_id"]),
                [
                    str(row["annotation_id"]),
                    int(row["taxid"]),
                    source_code(row["source_database"]),
                    matched,
                    _optional_text(row["release_date"]),
                    _optional_text(row["assembly_name"]),
                ],
            )
        )

    scored.sort(key=lambda t: (-t[0], t[1]))
    rows = [t[2] for t in scored]

    return _envelope(
        terms=terms_out, include=include, exclude=exclude, rows=rows
    )
