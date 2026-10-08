"""Cross-annotation hit queries over gene_corpus.sqlite."""
from __future__ import annotations

import heapq
import sqlite3
from collections.abc import Sequence

from helpers.curie import CurieError, format_curie, parse_curie
from helpers.cursor import (
    QueryError,
    decode_hits_annotations_cursor,
    decode_hits_cursor,
    encode_hits_annotations_cursor,
    encode_hits_cursor,
    filter_fingerprint,
    resolve_page_token,
)
from helpers.pagination import DEFAULT_LIMIT, MAX_CURIES_GET, clamp_limit
from helpers.strand import format_strand

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
        "strand": format_strand(row["strand"]),
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


def parse_curie_list(
    raw_values: str | Sequence[str],
    *,
    max_curies: int,
) -> list[str]:
    """Split comma-separated CURIEs, strip, and drop empty segments.

    Accepts one CSV string or an already-split sequence (each item may
    still contain commas). Requires at least one value and at most
    ``max_curies``.
    """
    pieces: Sequence[str] = (raw_values,) if isinstance(raw_values, str) else raw_values
    curies: list[str] = []
    for raw in pieces:
        if raw is None:
            continue
        for part in str(raw).split(","):
            text = part.strip()
            if text:
                curies.append(text)
    if not curies:
        raise QueryError(
            "at least one curie is required",
            code="curies_required",
        )
    if len(curies) > max_curies:
        raise QueryError(
            f"curies accepts at most {max_curies} values",
            code="too_many_curies",
            max=max_curies,
            received=len(curies),
        )
    return curies


def list_hits(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    taxid: int | None = None,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Any-match page of gene_hit rows for one or more CURIEs.

    No exact ``total``. Unknown prefixes become error entries; known
    prefixes with no rows contribute nothing.
    """
    curies = parse_curie_list(curies, max_curies=max_curies)
    if taxid is not None and int(taxid) < 1:
        raise QueryError(
            "taxid must be a positive integer",
            code="invalid_taxid",
        )

    limit = clamp_limit(limit)
    parsed: list[tuple[str, str, str]] = []  # (raw, namespace, accession)
    errors: list[dict] = []
    for raw in curies:
        try:
            namespace, accession = parse_curie(raw)
            parsed.append((raw, namespace, accession))
        except CurieError as exc:
            errors.append(
                {"curie": raw, "code": exc.code, "message": str(exc)}
            )

    # Deduplicate identical (namespace, accession) while keeping display ids.
    by_key: dict[tuple[str, str], list[str]] = {}
    for _raw, namespace, accession in parsed:
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )
    display_ids = sorted(
        {format_curie(ns, acc) for ns, acc in by_key}
    )
    filter_f = filter_fingerprint(
        {
            "curies": display_ids,
            "taxid": int(taxid) if taxid is not None else None,
        }
    )

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_ak: int | None = None
    bound_lid: int | None = None
    if page_token is not None:
        bound_ak, bound_lid = decode_hits_cursor(page_token, filter_f=filter_f)

    if taxid is not None and not taxid_has_annotations(conn, int(taxid)):
        return {
            "limit": limit,
            "results": [],
            "next": None,
            "previous": None,
            "errors": errors or None,
            "taxid": int(taxid),
        }

    after = None
    before = None
    if direction == "next" and bound_ak is not None and bound_lid is not None:
        after = (bound_ak, bound_lid)
    elif direction == "prev" and bound_ak is not None and bound_lid is not None:
        before = (bound_ak, bound_lid)

    # Heap merge: each stream contributes (ak, lid, display_ids, row).
    # Cap: at most MAX_CURIES_POST streams × (limit+1) rows — RAM-bounded.
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
                prev_token = encode_hits_cursor(
                    annotation_key=first_ak,
                    local_id=first_lid,
                    filter_f=filter_f,
                )
            next_token = encode_hits_cursor(
                annotation_key=last_ak,
                local_id=last_lid,
                filter_f=filter_f,
            )
        else:
            if has_more:
                next_token = encode_hits_cursor(
                    annotation_key=last_ak,
                    local_id=last_lid,
                    filter_f=filter_f,
                )
            if direction == "next":
                prev_token = encode_hits_cursor(
                    annotation_key=first_ak,
                    local_id=first_lid,
                    filter_f=filter_f,
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


def _fetch_annotation_groups(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
    *,
    after_ak: int | None,
    before_ak: int | None,
    limit: int,
    taxid: int | None = None,
) -> list[sqlite3.Row]:
    """Fetch up to *limit* annotation_key groups for one CURIE."""
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

    if after_ak is not None:
        seek = "AND h.annotation_key > ?"
        params.append(int(after_ak))
        order = "ORDER BY h.annotation_key ASC"
    elif before_ak is not None:
        seek = "AND h.annotation_key < ?"
        params.append(int(before_ak))
        order = "ORDER BY h.annotation_key DESC"
    else:
        seek = ""
        order = "ORDER BY h.annotation_key ASC"

    return conn.execute(
        f"""
        SELECT h.annotation_key AS annotation_key,
               COUNT(*) AS n_genes
        FROM gene_hit h
        WHERE h.namespace = ? AND h.accession = ?
        {lineage_sql}
        {seek}
        GROUP BY h.annotation_key
        {order}
        LIMIT ?
        """,
        [*params, limit],
    ).fetchall()


def list_hit_annotations(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    taxid: int | None = None,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Any-match page of annotations for one or more CURIEs (one row per annotation).

    No exact ``total``. Unknown prefixes become error entries.
    """
    curies = parse_curie_list(curies, max_curies=max_curies)
    if taxid is not None and int(taxid) < 1:
        raise QueryError(
            "taxid must be a positive integer",
            code="invalid_taxid",
        )

    limit = clamp_limit(limit)
    parsed: list[tuple[str, str, str]] = []
    errors: list[dict] = []
    for raw in curies:
        try:
            namespace, accession = parse_curie(raw)
            parsed.append((raw, namespace, accession))
        except CurieError as exc:
            errors.append(
                {"curie": raw, "code": exc.code, "message": str(exc)}
            )

    empty: dict = {
        "limit": limit,
        "results": [],
        "next": None,
        "previous": None,
    }
    if errors:
        empty["errors"] = errors
    if taxid is not None:
        empty["taxid"] = int(taxid)

    by_key: dict[tuple[str, str], list[str]] = {}
    for _raw, namespace, accession in parsed:
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )
    display_ids = sorted(
        {format_curie(ns, acc) for ns, acc in by_key}
    )
    filter_f = filter_fingerprint(
        {
            "curies": display_ids,
            "taxid": int(taxid) if taxid is not None else None,
        }
    )

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_ak: int | None = None
    if page_token is not None:
        bound_ak = decode_hits_annotations_cursor(
            page_token, filter_f=filter_f
        )

    if taxid is not None and not taxid_has_annotations(conn, int(taxid)):
        return empty

    after_ak = None
    before_ak = None
    if direction == "next" and bound_ak is not None:
        after_ak = bound_ak
    elif direction == "prev" and bound_ak is not None:
        before_ak = bound_ak

    taxid_i = int(taxid) if taxid is not None else None
    # Streams: (ak, display_ids)
    streams: list[list[tuple[int, list[str]]]] = []
    for (namespace, accession), display_ids_stream in by_key.items():
        rows = _fetch_annotation_groups(
            conn,
            namespace,
            accession,
            after_ak=after_ak,
            before_ak=before_ak,
            limit=limit + 1,
            taxid=taxid_i,
        )
        if before_ak is not None:
            rows = list(reversed(rows))
        stream = [
            (
                int(row["annotation_key"]),
                list(dict.fromkeys(display_ids_stream)),
            )
            for row in rows
        ]
        if stream:
            streams.append(stream)

    matched_by_ak: dict[int, list[str]] = {}
    order_keys: list[int] = []
    heap: list[tuple[int, int, int]] = []  # ak, stream_i, pos
    for i, stream in enumerate(streams):
        ak, _ = stream[0]
        heapq.heappush(heap, (ak, i, 0))

    while heap and len(order_keys) < limit + 1:
        ak, i, pos = heapq.heappop(heap)
        _, stream_display = streams[i][pos]
        if ak not in matched_by_ak:
            matched_by_ak[ak] = list(stream_display)
            order_keys.append(ak)
        else:
            seen = set(matched_by_ak[ak])
            for d in stream_display:
                if d not in seen:
                    matched_by_ak[ak].append(d)
                    seen.add(d)
        nxt = pos + 1
        if nxt < len(streams[i]):
            nak, _ = streams[i][nxt]
            heapq.heappush(heap, (nak, i, nxt))

    going_prev = direction == "prev"
    has_more = len(order_keys) > limit
    order_keys = order_keys[:limit]
    if going_prev:
        order_keys = list(reversed(order_keys))

    next_token = None
    prev_token = None
    if order_keys:
        first_ak = order_keys[0]
        last_ak = order_keys[-1]
        if going_prev:
            if has_more:
                prev_token = encode_hits_annotations_cursor(
                    annotation_key=first_ak,
                    filter_f=filter_f,
                )
            next_token = encode_hits_annotations_cursor(
                annotation_key=last_ak,
                filter_f=filter_f,
            )
        else:
            if has_more:
                next_token = encode_hits_annotations_cursor(
                    annotation_key=last_ak,
                    filter_f=filter_f,
                )
            if direction == "next":
                prev_token = encode_hits_annotations_cursor(
                    annotation_key=first_ak,
                    filter_f=filter_f,
                )

    if not order_keys:
        return empty

    # Exact n_genes under multi-CURIE overlap (page keys only).
    curie_preds = " OR ".join(
        "(h.namespace = ? AND h.accession = ?)" for _ in by_key
    )
    curie_params: list = []
    for namespace, accession in by_key:
        curie_params.extend([namespace, accession])
    ak_placeholders = ",".join("?" for _ in order_keys)
    lineage_sql = ""
    lineage_params: list = []
    if taxid_i is not None:
        lineage_sql = (
            "AND EXISTS ("
            "SELECT 1 FROM annotation_lineage al "
            "WHERE al.annotation_key = h.annotation_key AND al.taxid = ?"
            ")"
        )
        lineage_params.append(taxid_i)

    count_rows = conn.execute(
        f"""
        SELECT h.annotation_key AS annotation_key,
               COUNT(DISTINCT h.local_id) AS n_genes
        FROM gene_hit h
        WHERE h.annotation_key IN ({ak_placeholders})
          AND ({curie_preds})
          {lineage_sql}
        GROUP BY h.annotation_key
        """,
        [*order_keys, *curie_params, *lineage_params],
    ).fetchall()
    n_genes_by_ak = {
        int(row["annotation_key"]): int(row["n_genes"]) for row in count_rows
    }

    ann_rows = conn.execute(
        f"""
        SELECT annotation_key, annotation_id, taxid, assembly_accession,
               organism_name, source_database
        FROM annotation
        WHERE annotation_key IN ({ak_placeholders})
        """,
        order_keys,
    ).fetchall()
    ann_by_ak = {int(row["annotation_key"]): row for row in ann_rows}

    results = []
    for ak in order_keys:
        ann = ann_by_ak.get(ak)
        if ann is None:
            continue
        results.append(
            {
                "annotation_id": ann["annotation_id"],
                "taxid": int(ann["taxid"]),
                "assembly_accession": ann["assembly_accession"],
                "organism_name": ann["organism_name"],
                "source_database": ann["source_database"],
                "n_genes": n_genes_by_ak.get(ak, 0),
                "matched": matched_by_ak.get(ak, []),
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
