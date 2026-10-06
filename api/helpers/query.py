"""SQL queries over the locked global gene_corpus.sqlite serve schema."""
from __future__ import annotations

import heapq
import sqlite3
from collections.abc import Sequence

from helpers.curie import CurieError, format_curie, parse_curie
from helpers.cursor import (
    QueryError,
    decode_accessions_cursor,
    decode_genes_cursor,
    encode_accessions_cursor,
    encode_genes_cursor,
)

ACCESSION_SORTS = frozenset({"accession", "n_annotations", "n_loci"})
DEFAULT_LIMIT = 50
MAX_LIMIT = 100
MAX_IDS = 32


def _clamp_limit(limit: int) -> int:
    return max(1, min(int(limit), MAX_LIMIT))


def list_namespaces(conn: sqlite3.Connection) -> dict:
    rows = conn.execute(
        """
        SELECT namespace, accession_count
        FROM namespace_stats
        ORDER BY namespace
        """
    ).fetchall()
    results = [
        {
            "namespace": row["namespace"],
            "accession_count": int(row["accession_count"]),
        }
        for row in rows
    ]
    return {
        "total": len(results),
        "limit": len(results),
        "results": results,
        "next": None,
        "previous": None,
    }


def get_namespace(conn: sqlite3.Connection, namespace: str) -> dict | None:
    row = conn.execute(
        """
        SELECT namespace, accession_count
        FROM namespace_stats
        WHERE namespace = ?
        """,
        (namespace,),
    ).fetchone()
    if row is None:
        return None
    return {
        "namespace": row["namespace"],
        "accession_count": int(row["accession_count"]),
    }


def _accessions_page_sql(
    *,
    sort: str,
    sort_order: str,
    direction: str | None,
) -> tuple[str, str]:
    """Return (where_extra_sql, order_sql) for accession keyset pages."""
    descending = sort_order == "desc"
    count_col = "n_loci" if sort == "n_loci" else "n_annotations"

    if sort == "accession":
        if direction is None:
            order = "accession DESC" if descending else "accession ASC"
            return "", order
        if direction == "next":
            if descending:
                return "AND accession < ?", "accession DESC"
            return "AND accession > ?", "accession ASC"
        if descending:
            return "AND accession > ?", "accession ASC"
        return "AND accession < ?", "accession DESC"

    # sort == n_loci or n_annotations; stable tie-break accession ASC.
    if direction is None:
        if descending:
            return "", f"{count_col} DESC, accession ASC"
        return "", f"{count_col} ASC, accession ASC"

    if direction == "next":
        if descending:
            return (
                f"AND ({count_col} < ? OR ({count_col} = ? AND accession > ?))",
                f"{count_col} DESC, accession ASC",
            )
        return (
            f"AND ({count_col} > ? OR ({count_col} = ? AND accession > ?))",
            f"{count_col} ASC, accession ASC",
        )

    if descending:
        return (
            f"AND ({count_col} > ? OR ({count_col} = ? AND accession < ?))",
            f"{count_col} ASC, accession DESC",
        )
    return (
        f"AND ({count_col} < ? OR ({count_col} = ? AND accession < ?))",
        f"{count_col} DESC, accession DESC",
    )


def list_accessions(
    conn: sqlite3.Connection,
    namespace: str,
    *,
    sort: str = "n_loci",
    sort_order: str = "desc",
    cursor: str | None = None,
    limit: int = DEFAULT_LIMIT,
) -> dict | None:
    stats = get_namespace(conn, namespace)
    if stats is None:
        return None

    column = sort if sort in ACCESSION_SORTS else "n_loci"
    order = "asc" if (sort_order or "").lower() == "asc" else "desc"
    limit = _clamp_limit(limit)

    direction: str | None = None
    bound_n: int | None = None
    bound_a: str | None = None
    if cursor:
        direction, bound_n, bound_a = decode_accessions_cursor(
            cursor, sort=column, sort_order=order
        )

    where_extra, order_sql = _accessions_page_sql(
        sort=column, sort_order=order, direction=direction
    )
    params: list = [namespace]
    if direction is not None:
        assert bound_a is not None and bound_n is not None
        if column == "accession":
            params.append(bound_a)
        else:
            params.extend([bound_n, bound_n, bound_a])

    rows = conn.execute(
        f"""
        SELECT accession, n_annotations, n_loci
        FROM xref_meta
        WHERE namespace = ?
        {where_extra}
        ORDER BY {order_sql}
        LIMIT ?
        """,
        [*params, limit + 1],
    ).fetchall()

    going_prev = direction == "prev"
    has_more_in_scan = len(rows) > limit
    rows = rows[:limit]
    if going_prev:
        rows = list(reversed(rows))

    results = [
        {
            "accession": row["accession"],
            "n_annotations": int(row["n_annotations"]),
            "n_loci": int(row["n_loci"]),
        }
        for row in rows
    ]

    def _sort_value(item: dict) -> int:
        if column == "accession":
            return 0
        return int(item["n_loci" if column == "n_loci" else "n_annotations"])

    next_token = None
    prev_token = None
    if results:
        first = results[0]
        last = results[-1]
        if going_prev:
            if has_more_in_scan:
                prev_token = encode_accessions_cursor(
                    direction="prev",
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(first),
                    accession=first["accession"],
                )
            next_token = encode_accessions_cursor(
                direction="next",
                sort=column,
                sort_order=order,
                sort_value=_sort_value(last),
                accession=last["accession"],
            )
        else:
            if has_more_in_scan:
                next_token = encode_accessions_cursor(
                    direction="next",
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(last),
                    accession=last["accession"],
                )
            if direction == "next":
                prev_token = encode_accessions_cursor(
                    direction="prev",
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(first),
                    accession=first["accession"],
                )

    return {
        "namespace": namespace,
        "total": stats["accession_count"],
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }


def taxid_has_annotations(conn: sqlite3.Connection, taxid: int) -> bool:
    """True when any annotation lineage row carries *taxid*."""
    row = conn.execute(
        """
        SELECT 1 FROM annotation_lineage WHERE taxid = ? LIMIT 1
        """,
        (int(taxid),),
    ).fetchone()
    return row is not None


def _gene_row_dict(row: sqlite3.Row, *, matched: list[str]) -> dict:
    return {
        "annotation_id": row["annotation_id"],
        "taxid": int(row["taxid"]),
        "assembly_accession": row["assembly_accession"],
        "organism_name": row["organism_name"],
        "source_database": row["source_database"],
        "local_id": int(row["local_id"]),
        "seqid": row["seqid"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "strand": int(row["strand"]),
        "feature_type": row["feature_type"],
        "biotype": row["biotype"],
        "primary_name": row["primary_name"],
        "matched": matched,
        "_annotation_key": int(row["annotation_key"]),
    }


def _fetch_gene_page(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
    *,
    after: tuple[int, int] | None,
    before: tuple[int, int] | None,
    limit: int,
    taxid: int | None = None,
) -> list[sqlite3.Row]:
    """
    Fetch up to *limit* gene_hit rows for one CURIE.

    Optional *taxid* is applied in SQL via ``annotation_lineage`` so the
    engine skips non-matching keys without a Python key set or an unbounded
    batch-scan loop in RAM.
    """
    params: list = [namespace, accession]
    lineage_sql = ""
    if taxid is not None:
        lineage_sql = (
            "AND EXISTS ("
            "SELECT 1 FROM annotation_lineage al "
            "WHERE al.annotation_key = h.annotation_key AND al.taxid = ?"
            ")"
        )
        params.append(int(taxid))

    if after is not None:
        seek = (
            "AND (h.annotation_key > ? OR "
            "(h.annotation_key = ? AND h.local_id > ?))"
        )
        params.extend([after[0], after[0], after[1]])
        order = "ORDER BY h.annotation_key ASC, h.local_id ASC"
    elif before is not None:
        seek = (
            "AND (h.annotation_key < ? OR "
            "(h.annotation_key = ? AND h.local_id < ?))"
        )
        params.extend([before[0], before[0], before[1]])
        order = "ORDER BY h.annotation_key DESC, h.local_id DESC"
    else:
        seek = ""
        order = "ORDER BY h.annotation_key ASC, h.local_id ASC"

    return conn.execute(
        f"""
        SELECT
            h.annotation_key AS annotation_key,
            h.local_id AS local_id,
            h.seqid AS seqid,
            h.start AS start,
            h.end AS end,
            h.strand AS strand,
            h.feature_type AS feature_type,
            h.biotype AS biotype,
            h.primary_name AS primary_name,
            a.annotation_id AS annotation_id,
            a.taxid AS taxid,
            a.assembly_accession AS assembly_accession,
            a.organism_name AS organism_name,
            a.source_database AS source_database
        FROM gene_hit h
        JOIN annotation a ON a.annotation_key = h.annotation_key
        WHERE h.namespace = ? AND h.accession = ?
        {lineage_sql}
        {seek}
        {order}
        LIMIT ?
        """,
        [*params, limit],
    ).fetchall()


def list_genes(
    conn: sqlite3.Connection,
    ids: Sequence[str],
    *,
    taxid: int | None = None,
    cursor: str | None = None,
    limit: int = DEFAULT_LIMIT,
) -> dict:
    """
    Any-match page of gene_hit rows for one or more CURIEs.

    No exact ``total``. Unknown prefixes become error entries; known
    prefixes with no rows contribute nothing.
    """
    if not ids:
        raise QueryError("at least one id is required")
    if len(ids) > MAX_IDS:
        raise QueryError(f"id accepts at most {MAX_IDS} values")
    if taxid is not None and int(taxid) < 1:
        raise QueryError("taxid must be a positive integer")

    limit = _clamp_limit(limit)
    parsed: list[tuple[str, str, str]] = []  # (raw, namespace, accession)
    errors: list[dict] = []
    for raw in ids:
        try:
            namespace, accession = parse_curie(raw)
            parsed.append((raw.strip(), namespace, accession))
        except CurieError as exc:
            errors.append({"id": raw, "error": str(exc)})

    if taxid is not None and not taxid_has_annotations(conn, int(taxid)):
        return {
            "limit": limit,
            "results": [],
            "next": None,
            "previous": None,
            "errors": errors or None,
            "taxid": int(taxid),
        }

    direction: str | None = None
    bound_ak: int | None = None
    bound_lid: int | None = None
    if cursor:
        direction, bound_ak, bound_lid = decode_genes_cursor(cursor)

    # Deduplicate identical (namespace, accession) while keeping display ids.
    by_key: dict[tuple[str, str], list[str]] = {}
    for _raw, namespace, accession in parsed:
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )

    after = None
    before = None
    if direction == "next" and bound_ak is not None and bound_lid is not None:
        after = (bound_ak, bound_lid)
    elif direction == "prev" and bound_ak is not None and bound_lid is not None:
        before = (bound_ak, bound_lid)

    # Heap merge: each stream contributes (ak, lid, display_ids, row).
    # Cap: at most MAX_IDS streams × (limit+1) rows — RAM-bounded.
    streams: list[list[tuple[int, int, list[str], sqlite3.Row]]] = []
    taxid_i = int(taxid) if taxid is not None else None
    for (namespace, accession), display_ids in by_key.items():
        rows = _fetch_gene_page(
            conn,
            namespace,
            accession,
            after=after,
            before=before,
            limit=limit + 1,
            taxid=taxid_i,
        )
        if before is not None:
            rows = list(reversed(rows))
        stream = [
            (
                int(row["annotation_key"]),
                int(row["local_id"]),
                list(dict.fromkeys(display_ids)),
                row,
            )
            for row in rows
        ]
        if stream:
            streams.append(stream)

    merged: dict[tuple[int, int], dict] = {}
    order_keys: list[tuple[int, int]] = []
    heap: list[tuple[int, int, int, int]] = []
    for i, stream in enumerate(streams):
        ak, lid, _, _ = stream[0]
        heapq.heappush(heap, (ak, lid, i, 0))

    while heap and len(order_keys) < limit + 1:
        ak, lid, i, pos = heapq.heappop(heap)
        _, _, display_ids, row = streams[i][pos]
        key = (ak, lid)
        if key not in merged:
            merged[key] = _gene_row_dict(row, matched=list(display_ids))
            order_keys.append(key)
        else:
            seen = set(merged[key]["matched"])
            for d in display_ids:
                if d not in seen:
                    merged[key]["matched"].append(d)
                    seen.add(d)
        nxt = pos + 1
        if nxt < len(streams[i]):
            nak, nlid, _, _ = streams[i][nxt]
            heapq.heappush(heap, (nak, nlid, i, nxt))

    going_prev = direction == "prev"
    has_more = len(order_keys) > limit
    order_keys = order_keys[:limit]
    if going_prev:
        order_keys = list(reversed(order_keys))

    next_token = None
    prev_token = None
    if order_keys:
        first_ak, first_lid = order_keys[0]
        last_ak, last_lid = order_keys[-1]
        if going_prev:
            if has_more:
                prev_token = encode_genes_cursor(
                    direction="prev",
                    annotation_key=first_ak,
                    local_id=first_lid,
                )
            next_token = encode_genes_cursor(
                direction="next",
                annotation_key=last_ak,
                local_id=last_lid,
            )
        else:
            if has_more:
                next_token = encode_genes_cursor(
                    direction="next",
                    annotation_key=last_ak,
                    local_id=last_lid,
                )
            if direction == "next":
                prev_token = encode_genes_cursor(
                    direction="prev",
                    annotation_key=first_ak,
                    local_id=first_lid,
                )

    results = []
    for key in order_keys:
        item = merged[key]
        results.append(
            {
                "annotation_id": item["annotation_id"],
                "taxid": item["taxid"],
                "assembly_accession": item["assembly_accession"],
                "organism_name": item["organism_name"],
                "source_database": item["source_database"],
                "local_id": item["local_id"],
                "seqid": item["seqid"],
                "start": item["start"],
                "end": item["end"],
                "strand": item["strand"],
                "feature_type": item["feature_type"],
                "biotype": item["biotype"],
                "primary_name": item["primary_name"],
                "matched": item["matched"],
            }
        )

    out: dict = {
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }
    if errors:
        out["errors"] = errors
    if taxid is not None:
        out["taxid"] = int(taxid)
    return out
