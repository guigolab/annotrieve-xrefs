"""Cross-annotation hit queries over gene_corpus.sqlite."""
from __future__ import annotations

import heapq
import sqlite3
from collections.abc import Sequence
from typing import Literal

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

MATCH_ANY = "any"
MATCH_ALL = "all"
MATCH_MODES = frozenset({MATCH_ANY, MATCH_ALL})

# Batch size when scanning the driver CURIE under match=all (many rows fail probes).
_ALL_DRIVER_BATCH = 256

MatchMode = Literal["any", "all"]


def parse_match(raw: str | None) -> MatchMode:
    """Validate ``match``; default ``any``."""
    text = MATCH_ANY if raw is None else str(raw).strip().lower()
    if text not in MATCH_MODES:
        raise QueryError(
            "match must be any or all",
            code="invalid_match",
        )
    return text  # type: ignore[return-value]


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


def _lineage_sql(taxid: int | None, *, alias: str = "h") -> tuple[str, list]:
    if taxid is None:
        return "", []
    return (
        (
            f"AND EXISTS ("
            f"SELECT 1 FROM annotation_lineage al "
            f"WHERE al.annotation_key = {alias}.annotation_key AND al.taxid = ?"
            f")"
        ),
        [int(taxid)],
    )


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
    lineage_sql, lineage_params = _lineage_sql(taxid)
    params.extend(lineage_params)

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


def _parse_hits_request(
    curies: str | Sequence[str],
    *,
    taxid: int | None,
    match: str | None,
    max_curies: int,
) -> tuple[
    MatchMode,
    list[dict],
    dict[tuple[str, str], list[str]],
    list[str],
    str,
    int | None,
]:
    """
    Shared request parsing for hits routes.

    Returns ``(match, errors, by_key, display_ids, filter_f, taxid_i)``.
    """
    match_mode = parse_match(match)
    curies = parse_curie_list(curies, max_curies=max_curies)
    if taxid is not None and int(taxid) < 1:
        raise QueryError(
            "taxid must be a positive integer",
            code="invalid_taxid",
        )

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

    by_key: dict[tuple[str, str], list[str]] = {}
    for _raw, namespace, accession in parsed:
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )
    display_ids = sorted(
        {format_curie(ns, acc) for ns, acc in by_key}
    )
    taxid_i = int(taxid) if taxid is not None else None
    filter_f = filter_fingerprint(
        {
            "curies": display_ids,
            "taxid": taxid_i,
            "match": match_mode,
        }
    )
    return match_mode, errors, by_key, display_ids, filter_f, taxid_i


def _xref_meta_n_loci(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
) -> int | None:
    """Return ``n_loci`` from xref_meta, or None when the accession is absent."""
    row = conn.execute(
        """
        SELECT n_loci FROM xref_meta
        WHERE namespace = ? AND accession = ?
        """,
        (namespace, accession),
    ).fetchone()
    if row is None:
        return None
    return int(row["n_loci"])


def _pick_driver(
    conn: sqlite3.Connection,
    by_key: dict[tuple[str, str], list[str]],
) -> tuple[str, str] | None:
    """
    Smallest CURIE by ``xref_meta.n_loci``; ties by ``(namespace, accession)``.

    Returns None when any key is missing from xref_meta (empty intersection).
    """
    best: tuple[int, str, str] | None = None
    for namespace, accession in by_key:
        n_loci = _xref_meta_n_loci(conn, namespace, accession)
        if n_loci is None:
            return None
        cand = (n_loci, namespace, accession)
        if best is None or cand < best:
            best = cand
    if best is None:
        return None
    return best[1], best[2]


def _probe_gene_hit(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
    annotation_key: int,
    local_id: int,
) -> bool:
    row = conn.execute(
        """
        SELECT 1 FROM gene_hit
        WHERE namespace = ? AND accession = ?
          AND annotation_key = ? AND local_id = ?
        LIMIT 1
        """,
        (namespace, accession, int(annotation_key), int(local_id)),
    ).fetchone()
    return row is not None


def _probe_annotation_hit(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
    annotation_key: int,
) -> bool:
    row = conn.execute(
        """
        SELECT 1 FROM gene_hit
        WHERE namespace = ? AND accession = ? AND annotation_key = ?
        LIMIT 1
        """,
        (namespace, accession, int(annotation_key)),
    ).fetchone()
    return row is not None


def _all_matched_display_ids(
    by_key: dict[tuple[str, str], list[str]],
) -> list[str]:
    """Stable unique display CURIEs for a full intersection hit."""
    out: list[str] = []
    seen: set[str] = set()
    for display_list in by_key.values():
        for d in display_list:
            if d not in seen:
                seen.add(d)
                out.append(d)
    return sorted(out)


def _hits_empty(
    *,
    limit: int,
    errors: list[dict],
    taxid: int | None,
) -> dict:
    out: dict = {
        "limit": limit,
        "results": [],
        "next": None,
        "previous": None,
    }
    if errors:
        out["errors"] = errors
    if taxid is not None:
        out["taxid"] = int(taxid)
    return out


def _page_tokens_hits(
    order_keys: list[tuple[int, int]],
    *,
    direction: str | None,
    has_more: bool,
    filter_f: str,
) -> tuple[str | None, str | None]:
    next_token = None
    prev_token = None
    if not order_keys:
        return None, None
    first_ak, first_lid = order_keys[0]
    last_ak, last_lid = order_keys[-1]
    going_prev = direction == "prev"
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
    return next_token, prev_token


def _page_tokens_annotations(
    order_keys: list[int],
    *,
    direction: str | None,
    has_more: bool,
    filter_f: str,
) -> tuple[str | None, str | None]:
    next_token = None
    prev_token = None
    if not order_keys:
        return None, None
    first_ak = order_keys[0]
    last_ak = order_keys[-1]
    going_prev = direction == "prev"
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
    return next_token, prev_token


def _list_hits_any(
    conn: sqlite3.Connection,
    *,
    by_key: dict[tuple[str, str], list[str]],
    errors: list[dict],
    filter_f: str,
    taxid_i: int | None,
    direction: str | None,
    after: tuple[int, int] | None,
    before: tuple[int, int] | None,
    limit: int,
) -> dict:
    streams: list[list[tuple[int, int, list[str], sqlite3.Row]]] = []
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

    next_token, prev_token = _page_tokens_hits(
        order_keys,
        direction=direction,
        has_more=has_more,
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
    if taxid_i is not None:
        out["taxid"] = taxid_i
    return out


def _list_hits_all(
    conn: sqlite3.Connection,
    *,
    by_key: dict[tuple[str, str], list[str]],
    errors: list[dict],
    filter_f: str,
    taxid_i: int | None,
    direction: str | None,
    after: tuple[int, int] | None,
    before: tuple[int, int] | None,
    limit: int,
) -> dict:
    """Intersection: scan smallest CURIE, probe others on the primary key."""
    empty = _hits_empty(limit=limit, errors=errors, taxid=taxid_i)
    if errors or not by_key:
        return empty

    driver = _pick_driver(conn, by_key)
    if driver is None:
        return empty

    driver_ns, driver_acc = driver
    others = [(ns, acc) for (ns, acc) in by_key if (ns, acc) != driver]
    matched_all = _all_matched_display_ids(by_key)

    going_prev = direction == "prev"
    collected: list[sqlite3.Row] = []
    # Seek past the last accepted key while paging forward/back.
    seek_after = after
    seek_before = before

    while len(collected) < limit + 1:
        if going_prev:
            batch = _fetch_gene_page(
                conn,
                driver_ns,
                driver_acc,
                after=None,
                before=seek_before,
                limit=_ALL_DRIVER_BATCH,
                taxid=taxid_i,
            )
        else:
            batch = _fetch_gene_page(
                conn,
                driver_ns,
                driver_acc,
                after=seek_after,
                before=None,
                limit=_ALL_DRIVER_BATCH,
                taxid=taxid_i,
            )
        if not batch:
            break

        for row in batch:
            ak = int(row["annotation_key"])
            lid = int(row["local_id"])
            if all(
                _probe_gene_hit(conn, ns, acc, ak, lid) for ns, acc in others
            ):
                collected.append(row)
                if len(collected) >= limit + 1:
                    break

        last = batch[-1]
        if going_prev:
            seek_before = (int(last["annotation_key"]), int(last["local_id"]))
        else:
            seek_after = (int(last["annotation_key"]), int(last["local_id"]))
        if len(batch) < _ALL_DRIVER_BATCH:
            break

    has_more = len(collected) > limit
    page_rows = collected[:limit]
    if going_prev:
        page_rows = list(reversed(page_rows))

    order_keys = [
        (int(row["annotation_key"]), int(row["local_id"])) for row in page_rows
    ]
    next_token, prev_token = _page_tokens_hits(
        order_keys,
        direction=direction,
        has_more=has_more,
        filter_f=filter_f,
    )

    results = []
    for row in page_rows:
        item = _gene_row_dict(row, matched=list(matched_all))
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
    if taxid_i is not None:
        out["taxid"] = taxid_i
    return out


def list_hits(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    taxid: int | None = None,
    match: str | None = MATCH_ANY,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Page of gene_hit rows for one or more CURIEs.

    ``match=any`` (default): union via heap merge. ``match=all``: intersection
    via smallest-CURIE scan and primary-key probes.

    No exact ``total``. Unknown prefixes become error entries; under ``all``
    any error or missing xref_meta accession yields an empty page.
    """
    match_mode, errors, by_key, _display_ids, filter_f, taxid_i = (
        _parse_hits_request(
            curies, taxid=taxid, match=match, max_curies=max_curies
        )
    )
    limit = clamp_limit(limit)

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_ak: int | None = None
    bound_lid: int | None = None
    if page_token is not None:
        bound_ak, bound_lid = decode_hits_cursor(page_token, filter_f=filter_f)

    if taxid_i is not None and not taxid_has_annotations(conn, taxid_i):
        return _hits_empty(limit=limit, errors=errors, taxid=taxid_i)

    after = None
    before = None
    if direction == "next" and bound_ak is not None and bound_lid is not None:
        after = (bound_ak, bound_lid)
    elif direction == "prev" and bound_ak is not None and bound_lid is not None:
        before = (bound_ak, bound_lid)

    # Single unique CURIE: intersection equals union.
    if match_mode == MATCH_ALL and (errors or len(by_key) > 1):
        return _list_hits_all(
            conn,
            by_key=by_key,
            errors=errors,
            filter_f=filter_f,
            taxid_i=taxid_i,
            direction=direction,
            after=after,
            before=before,
            limit=limit,
        )

    return _list_hits_any(
        conn,
        by_key=by_key,
        errors=errors,
        filter_f=filter_f,
        taxid_i=taxid_i,
        direction=direction,
        after=after,
        before=before,
        limit=limit,
    )


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
    lineage_sql, lineage_params = _lineage_sql(taxid)
    params.extend(lineage_params)

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


def _next_driver_annotation_key(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
    *,
    after_ak: int | None,
    before_ak: int | None,
    taxid: int | None,
) -> int | None:
    """Seek one annotation_key on the driver CURIE (ascending or descending)."""
    params: list = [namespace, accession]
    lineage_sql, lineage_params = _lineage_sql(taxid)
    params.extend(lineage_params)

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

    row = conn.execute(
        f"""
        SELECT h.annotation_key AS annotation_key
        FROM gene_hit h
        WHERE h.namespace = ? AND h.accession = ?
        {lineage_sql}
        {seek}
        {order}
        LIMIT 1
        """,
        params,
    ).fetchone()
    if row is None:
        return None
    return int(row["annotation_key"])


def _build_annotation_results(
    conn: sqlite3.Connection,
    *,
    order_keys: list[int],
    by_key: dict[tuple[str, str], list[str]],
    matched_by_ak: dict[int, list[str]],
    taxid_i: int | None,
) -> list[dict]:
    if not order_keys:
        return []

    curie_preds = " OR ".join(
        "(h.namespace = ? AND h.accession = ?)" for _ in by_key
    )
    curie_params: list = []
    for namespace, accession in by_key:
        curie_params.extend([namespace, accession])
    ak_placeholders = ",".join("?" for _ in order_keys)
    lineage_sql, lineage_params = _lineage_sql(taxid_i)

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
    return results


def _list_hit_annotations_any(
    conn: sqlite3.Connection,
    *,
    by_key: dict[tuple[str, str], list[str]],
    errors: list[dict],
    filter_f: str,
    taxid_i: int | None,
    direction: str | None,
    after_ak: int | None,
    before_ak: int | None,
    limit: int,
) -> dict:
    empty = _hits_empty(limit=limit, errors=errors, taxid=taxid_i)

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
    heap: list[tuple[int, int, int]] = []
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

    next_token, prev_token = _page_tokens_annotations(
        order_keys,
        direction=direction,
        has_more=has_more,
        filter_f=filter_f,
    )

    if not order_keys:
        return empty

    results = _build_annotation_results(
        conn,
        order_keys=order_keys,
        by_key=by_key,
        matched_by_ak=matched_by_ak,
        taxid_i=taxid_i,
    )

    out: dict = {
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }
    if errors:
        out["errors"] = errors
    if taxid_i is not None:
        out["taxid"] = taxid_i
    return out


def _list_hit_annotations_all(
    conn: sqlite3.Connection,
    *,
    by_key: dict[tuple[str, str], list[str]],
    errors: list[dict],
    filter_f: str,
    taxid_i: int | None,
    direction: str | None,
    after_ak: int | None,
    before_ak: int | None,
    limit: int,
) -> dict:
    """Intersection at annotation grain: seek driver keys, probe other CURIEs."""
    empty = _hits_empty(limit=limit, errors=errors, taxid=taxid_i)
    if errors or not by_key:
        return empty

    driver = _pick_driver(conn, by_key)
    if driver is None:
        return empty

    driver_ns, driver_acc = driver
    others = [(ns, acc) for (ns, acc) in by_key if (ns, acc) != driver]
    matched_all = _all_matched_display_ids(by_key)

    going_prev = direction == "prev"
    seek_after = after_ak
    seek_before = before_ak
    collected: list[int] = []

    while len(collected) < limit + 1:
        if going_prev:
            ak = _next_driver_annotation_key(
                conn,
                driver_ns,
                driver_acc,
                after_ak=None,
                before_ak=seek_before,
                taxid=taxid_i,
            )
        else:
            ak = _next_driver_annotation_key(
                conn,
                driver_ns,
                driver_acc,
                after_ak=seek_after,
                before_ak=None,
                taxid=taxid_i,
            )
        if ak is None:
            break
        if all(_probe_annotation_hit(conn, ns, acc, ak) for ns, acc in others):
            collected.append(ak)
        if going_prev:
            seek_before = ak
        else:
            seek_after = ak

    has_more = len(collected) > limit
    order_keys = collected[:limit]
    if going_prev:
        order_keys = list(reversed(order_keys))

    next_token, prev_token = _page_tokens_annotations(
        order_keys,
        direction=direction,
        has_more=has_more,
        filter_f=filter_f,
    )

    if not order_keys:
        return empty

    matched_by_ak = {ak: list(matched_all) for ak in order_keys}
    results = _build_annotation_results(
        conn,
        order_keys=order_keys,
        by_key=by_key,
        matched_by_ak=matched_by_ak,
        taxid_i=taxid_i,
    )

    out: dict = {
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }
    if errors:
        out["errors"] = errors
    if taxid_i is not None:
        out["taxid"] = taxid_i
    return out


def list_hit_annotations(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    taxid: int | None = None,
    match: str | None = MATCH_ANY,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Page of annotations for one or more CURIEs (one row per annotation).

    ``match=any`` (default): union. ``match=all``: every CURIE must appear in
    the annotation (on any of its genes). ``n_genes`` is the distinct count of
    genes that matched at least one CURIE on the page.
    """
    match_mode, errors, by_key, _display_ids, filter_f, taxid_i = (
        _parse_hits_request(
            curies, taxid=taxid, match=match, max_curies=max_curies
        )
    )
    limit = clamp_limit(limit)
    empty = _hits_empty(limit=limit, errors=errors, taxid=taxid_i)

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_ak: int | None = None
    if page_token is not None:
        bound_ak = decode_hits_annotations_cursor(
            page_token, filter_f=filter_f
        )

    if taxid_i is not None and not taxid_has_annotations(conn, taxid_i):
        return empty

    after_ak = None
    before_ak = None
    if direction == "next" and bound_ak is not None:
        after_ak = bound_ak
    elif direction == "prev" and bound_ak is not None:
        before_ak = bound_ak

    if match_mode == MATCH_ALL and (errors or len(by_key) > 1):
        return _list_hit_annotations_all(
            conn,
            by_key=by_key,
            errors=errors,
            filter_f=filter_f,
            taxid_i=taxid_i,
            direction=direction,
            after_ak=after_ak,
            before_ak=before_ak,
            limit=limit,
        )

    return _list_hit_annotations_any(
        conn,
        by_key=by_key,
        errors=errors,
        filter_f=filter_f,
        taxid_i=taxid_i,
        direction=direction,
        after_ak=after_ak,
        before_ak=before_ak,
        limit=limit,
    )
