"""Annotation shard gene / xref queries (per-GFF genes.sqlite)."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from pathlib import Path

from helpers.curie import (
    CASEFOLD_NAMESPACES,
    CurieError,
    parse_curie_local,
    resolve_namespace_local,
)
from helpers.cursor import QueryError
from helpers.match import (
    MATCH_ALL,
    MATCH_ANY,
    MatchMode,
    parse_curie_list,
    parse_match,
)
from helpers.names import parse_exact_name, parse_name_list
from helpers.paths import per_gff_root, shard_genes_path
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_CURIES_GET,
    MAX_NAMES_GET,
    MIN_Q_LEN,
    clamp_limit,
    clamp_offset,
    page_envelope,
)
from helpers.session import GeneCorpusUnavailable
from helpers.strand import format_strand

_SHARD_CACHE_SIZE = -16384  # 16 MiB

_GENE_COLS = (
    "local_id, source_gene_id, feature_type, seqid, start, end, "
    "strand, biotype, primary_name, prose, prose_kind"
)


class ShardUnavailable(Exception):
    """Raised when a per-GFF genes.sqlite is missing."""


def open_shard(
    annotation_id: str,
    *,
    base: Path | None = None,
) -> sqlite3.Connection:
    """
    Open one annotation's genes.sqlite read-only.

    Caller must close the connection. Not process-cached.
    """
    if per_gff_root(base) is None:
        raise GeneCorpusUnavailable("Gene corpus is not published")
    path = shard_genes_path(annotation_id, base=base)
    if path is None or not path.is_file():
        raise ShardUnavailable(
            f"Gene shard not found for annotation {annotation_id}"
        )
    uri = f"file:{path.resolve()}?mode=ro&immutable=1"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA cache_size={_SHARD_CACHE_SIZE}")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA mmap_size=0")
    conn.execute("PRAGMA query_only=ON")
    return conn


def _gene_card(row: sqlite3.Row) -> dict:
    return {
        "local_id": int(row["local_id"]),
        "source_gene_id": row["source_gene_id"],
        "feature_type": row["feature_type"],
        "seqid": row["seqid"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "strand": format_strand(row["strand"]),
        "biotype": row["biotype"],
        "primary_name": row["primary_name"],
        "prose": row["prose"],
        "prose_kind": row["prose_kind"],
    }


def get_gene(
    conn: sqlite3.Connection,
    local_id: int,
) -> dict | None:
    """Return one gene card, or None if missing."""
    if local_id < 0:
        return None
    row = conn.execute(
        f"""
        SELECT {_GENE_COLS}
        FROM gene
        WHERE local_id = ?
        """,
        (int(local_id),),
    ).fetchone()
    if row is None:
        return None
    return _gene_card(row)


def gene_exists(conn: sqlite3.Connection, local_id: int) -> bool:
    if local_id < 0:
        return False
    row = conn.execute(
        "SELECT 1 FROM gene WHERE local_id = ? LIMIT 1",
        (int(local_id),),
    ).fetchone()
    return row is not None


def list_gene_xrefs(
    conn: sqlite3.Connection,
    annotation_id: str,
    local_id: int,
    *,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
) -> dict | None:
    """Paginated xref chips for one gene. None when the gene is missing."""
    if not gene_exists(conn, local_id):
        return None

    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    total = int(
        conn.execute(
            "SELECT COUNT(*) FROM gene_xref WHERE local_id = ?",
            (int(local_id),),
        ).fetchone()[0]
    )
    rows = conn.execute(
        """
        SELECT namespace, accession, display, origin, via_level
        FROM gene_xref
        WHERE local_id = ?
        ORDER BY namespace ASC, accession ASC
        LIMIT ? OFFSET ?
        """,
        (int(local_id), limit, offset),
    ).fetchall()
    results = [
        {
            "namespace": row["namespace"],
            "accession": row["accession"],
            "display": row["display"],
            "origin": row["origin"],
            "via_level": int(row["via_level"]),
        }
        for row in rows
    ]
    # annotation_id / local_id are path context only; not echoed.
    _ = annotation_id
    return page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
    )


def list_annotation_namespaces(
    conn: sqlite3.Connection,
    annotation_id: str,
) -> dict:
    """
    All namespaces present on one annotation shard (from gene_xref).

    ``accession_count`` is the number of distinct accessions per namespace.
    """
    # annotation_id is path context only; not echoed in the payload.
    _ = annotation_id
    rows = conn.execute(
        """
        SELECT namespace, COUNT(DISTINCT accession) AS accession_count
        FROM gene_xref
        GROUP BY namespace
        ORDER BY namespace
        """
    ).fetchall()
    results = [
        {
            "namespace": str(row["namespace"]),
            "accession_count": int(row["accession_count"]),
        }
        for row in rows
    ]
    return page_envelope(
        total=len(results),
        limit=len(results),
        offset=0,
        results=results,
    )


def resolve_path_namespace(raw: str) -> str | None:
    """Map a path ``{namespace}`` segment to the stored namespace form."""
    return resolve_namespace_local(raw)


def namespace_present(conn: sqlite3.Connection, namespace: str) -> bool:
    """True when ``namespace`` has at least one gene_xref row on this shard."""
    row = conn.execute(
        "SELECT 1 FROM gene_xref WHERE namespace = ? LIMIT 1",
        (namespace,),
    ).fetchone()
    return row is not None


def list_namespace_accessions(
    conn: sqlite3.Connection,
    annotation_id: str,
    namespace: str,
    *,
    q: str | None = None,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
) -> dict | None:
    """
    Paginated accessions under one namespace on this shard.

    Returns None when the namespace is absent. Each result:
    ``{accession, display, gene_count}``.
    """
    _ = annotation_id
    if not namespace_present(conn, namespace):
        return None

    limit = clamp_limit(limit)
    offset = clamp_offset(offset)

    prefix = (q or "").strip()
    if prefix and namespace in CASEFOLD_NAMESPACES:
        prefix = prefix.casefold()
    params: list = [namespace]
    where = "namespace = ?"
    if prefix:
        if len(prefix) < MIN_Q_LEN:
            raise QueryError(
                f"q requires at least {MIN_Q_LEN} characters",
                code="prefix_too_short",
                min_length=MIN_Q_LEN,
            )
        where += " AND accession >= ? AND accession < ?"
        params.extend([prefix, prefix + "\uffff"])

    total = int(
        conn.execute(
            f"""
            SELECT COUNT(DISTINCT accession)
            FROM gene_xref
            WHERE {where}
            """,
            params,
        ).fetchone()[0]
    )
    rows = conn.execute(
        f"""
        SELECT accession,
               MAX(display) AS display,
               COUNT(DISTINCT local_id) AS gene_count
        FROM gene_xref
        WHERE {where}
        GROUP BY accession
        ORDER BY accession ASC
        LIMIT ? OFFSET ?
        """,
        [*params, limit, offset],
    ).fetchall()
    results = [
        {
            "accession": str(row["accession"]),
            "display": row["display"],
            "gene_count": int(row["gene_count"]),
        }
        for row in rows
    ]
    return page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
    )


def _local_ids_for_accessions(
    conn: sqlite3.Connection, accessions: list[str]
) -> set[int]:
    """local_ids with any symbol/alias accession in ``accessions``."""
    if not accessions:
        return set()
    placeholders = ",".join("?" for _ in accessions)
    rows = conn.execute(
        f"""
        SELECT DISTINCT local_id
        FROM gene_xref
        WHERE namespace IN ('symbol', 'alias')
          AND accession IN ({placeholders})
        """,
        accessions,
    ).fetchall()
    return {int(row["local_id"]) for row in rows}


def _ids_by_keys(
    conn: sqlite3.Connection,
    keys: list[tuple[str, str]],
) -> dict[tuple[str, str], set[int]]:
    """
    Map each exact ``(namespace, accession)`` to its local_id set.

    One round-trip using ``xref_resolve``-friendly equality OR clauses.
    """
    out: dict[tuple[str, str], set[int]] = {key: set() for key in keys}
    if not keys:
        return out
    unique = list(dict.fromkeys(keys))
    clauses: list[str] = []
    params: list = []
    for ns, acc in unique:
        clauses.append("(namespace = ? AND accession = ?)")
        params.extend([ns, acc])
    rows = conn.execute(
        f"""
        SELECT namespace, accession, local_id
        FROM gene_xref
        WHERE {" OR ".join(clauses)}
        """,
        params,
    ).fetchall()
    for row in rows:
        key = (str(row["namespace"]), str(row["accession"]))
        if key in out:
            out[key].add(int(row["local_id"]))
    return out


def _combine_term_ids(
    per_term: list[set[int]],
    match: MatchMode,
) -> set[int]:
    if not per_term:
        return set()
    if match == MATCH_ALL:
        local_ids = set(per_term[0])
        for extra in per_term[1:]:
            local_ids &= extra
        return local_ids
    local_ids: set[int] = set()
    for ids in per_term:
        local_ids |= ids
    return local_ids


def _resolve_name_terms(
    conn: sqlite3.Connection,
    names: str | Sequence[str],
    *,
    match: MatchMode,
    max_names: int,
) -> tuple[list[dict], set[int]]:
    """
    Resolve exact ``names`` into term specs and the filtered local_id set.

    Each term spec: ``{input, candidates}`` (one casefolded accession).
    Wildcards are rejected (``wildcard_not_supported``).
    """
    raw_names = parse_name_list(names, max_names=max_names)
    terms: list[dict] = []
    per_term: list[set[int]] = []
    for raw in raw_names:
        inp, name_form = parse_exact_name(raw)
        candidates = [name_form]
        terms.append({"input": inp, "candidates": candidates})
        per_term.append(_local_ids_for_accessions(conn, candidates))

    return terms, _combine_term_ids(per_term, match)


def _resolve_curie_terms(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    match: MatchMode,
    max_curies: int,
) -> tuple[list[dict], set[int], list[set[int]]]:
    """
    Parse ``curies`` into term specs, combined local_ids, and per-term sets.

    Each term: ``{input, ns, accession}``. Exact keys only. One SQL round-trip.
    """
    raw_curies = parse_curie_list(curies, max_curies=max_curies)
    terms: list[dict] = []
    for raw in raw_curies:
        try:
            ns, acc = parse_curie_local(raw)
        except CurieError as exc:
            raise QueryError(str(exc), code=exc.code) from exc
        terms.append({"input": raw, "ns": ns, "accession": acc})

    keys = [(t["ns"], t["accession"]) for t in terms]
    by_key = _ids_by_keys(conn, keys)
    per_term = [by_key[(t["ns"], t["accession"])] for t in terms]
    return terms, _combine_term_ids(per_term, match), per_term


def _matches_for_page(
    conn: sqlite3.Connection,
    page_ids: list[int],
    terms: list[dict],
) -> dict[int, list[dict]]:
    """
    Map page local_id → inline matches ``{input, name, ns}``.

    Order per gene: term index, then namespace, then accession.
    """
    out: dict[int, list[dict]] = {lid: [] for lid in page_ids}
    if not page_ids or not terms:
        return out

    all_candidates = list(
        dict.fromkeys(
            acc for term in terms for acc in term["candidates"]
        )
    )
    if not all_candidates:
        return out

    # accession → terms that claimed it (term-list order)
    acc_to_terms: dict[str, list[dict]] = {}
    for term in terms:
        for acc in term["candidates"]:
            acc_to_terms.setdefault(acc, []).append(term)

    id_ph = ",".join("?" for _ in page_ids)
    acc_ph = ",".join("?" for _ in all_candidates)
    rows = conn.execute(
        f"""
        SELECT local_id, namespace, accession
        FROM gene_xref
        WHERE local_id IN ({id_ph})
          AND namespace IN ('symbol', 'alias')
          AND accession IN ({acc_ph})
        """,
        [*page_ids, *all_candidates],
    ).fetchall()

    # Collect then sort for stable output
    raw: dict[int, list[tuple[int, str, str, str]]] = {
        lid: [] for lid in page_ids
    }
    term_index = {id(term): i for i, term in enumerate(terms)}
    for row in rows:
        lid = int(row["local_id"])
        if lid not in raw:
            continue
        acc = str(row["accession"])
        ns = str(row["namespace"])
        for term in acc_to_terms.get(acc, []):
            raw[lid].append(
                (term_index[id(term)], ns, acc, term["input"])
            )

    for lid, items in raw.items():
        items.sort(key=lambda t: (t[0], t[1], t[2]))
        # Dedupe identical (term, ns, name) while keeping order
        seen: set[tuple[int, str, str]] = set()
        matches: list[dict] = []
        for ti, ns, name, inp in items:
            key = (ti, ns, name)
            if key in seen:
                continue
            seen.add(key)
            matches.append({"input": inp, "name": name, "ns": ns})
        out[lid] = matches
    return out


def _curie_matches_for_page(
    page_ids: list[int],
    terms: list[dict],
    per_term: list[set[int]],
) -> dict[int, list[dict]]:
    """
    Map page local_id → inline matches ``{input, ns, accession}``.

    Built from resolve-time per-term id sets (no extra SQL). Order per gene
    follows term list order; duplicates of the same (term, ns, accession)
    are dropped.
    """
    out: dict[int, list[dict]] = {lid: [] for lid in page_ids}
    if not page_ids or not terms:
        return out
    page_set = set(page_ids)
    for term, ids in zip(terms, per_term, strict=True):
        match = {
            "input": term["input"],
            "ns": term["ns"],
            "accession": term["accession"],
        }
        for lid in ids:
            if lid in page_set:
                out[lid].append(match)
    return out


def list_annotation_genes(
    conn: sqlite3.Connection,
    annotation_id: str,
    *,
    names: str | Sequence[str] | None = None,
    curies: str | Sequence[str] | None = None,
    match: str | None = MATCH_ANY,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
    max_names: int = MAX_NAMES_GET,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Paginated gene cards for one annotation shard.

    Optional ``names``: exact symbol∪alias filter (``match=any`` / ``match=all``).
    Wildcards are rejected (``wildcard_not_supported``).
    Optional ``curies``: exact ``ns:accession`` keys (any shard namespace).
    ``names`` and ``curies`` are mutually exclusive.
    When either filter is set, each gene card includes ``matches``.
    """
    # annotation_id is path context only; not echoed in the payload.
    _ = annotation_id
    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    match_mode = parse_match(match)

    if names is not None and curies is not None:
        raise QueryError(
            "names and curies cannot be combined",
            code="names_and_curies",
        )

    empty = page_envelope(
        total=0,
        limit=limit,
        offset=offset,
        results=[],
    )

    local_ids: set[int] | None = None
    name_terms: list[dict] = []
    curie_terms: list[dict] = []
    curie_per_term: list[set[int]] = []
    if names is not None:
        name_terms, local_ids = _resolve_name_terms(
            conn, names, match=match_mode, max_names=max_names
        )
        if not local_ids:
            return empty
    elif curies is not None:
        curie_terms, local_ids, curie_per_term = _resolve_curie_terms(
            conn, curies, match=match_mode, max_curies=max_curies
        )
        if not local_ids:
            return empty

    params: list = []
    where_parts: list[str] = ["1=1"]

    if local_ids is not None:
        placeholders = ",".join("?" for _ in local_ids)
        params.extend(sorted(local_ids))
        where_parts.append(f"g.local_id IN ({placeholders})")

    where_sql = " AND ".join(where_parts)

    total = int(
        conn.execute(
            f"""
            SELECT COUNT(*) FROM gene g
            WHERE {where_sql}
            """,
            params,
        ).fetchone()[0]
    )

    rows = conn.execute(
        f"""
        SELECT g.local_id, g.source_gene_id, g.feature_type, g.seqid,
               g.start, g.end, g.strand, g.biotype, g.primary_name,
               g.prose, g.prose_kind
        FROM gene g
        WHERE {where_sql}
        ORDER BY g.local_id ASC
        LIMIT ? OFFSET ?
        """,
        [*params, limit, offset],
    ).fetchall()

    results = [_gene_card(row) for row in rows]
    if results and name_terms:
        page_ids = [int(card["local_id"]) for card in results]
        by_lid = _matches_for_page(conn, page_ids, name_terms)
        for card in results:
            card["matches"] = by_lid.get(int(card["local_id"]), [])
    elif results and curie_terms:
        page_ids = [int(card["local_id"]) for card in results]
        by_lid = _curie_matches_for_page(
            page_ids, curie_terms, curie_per_term
        )
        for card in results:
            card["matches"] = by_lid.get(int(card["local_id"]), [])

    return page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
    )
