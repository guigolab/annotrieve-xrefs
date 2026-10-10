"""Annotation shard gene / xref queries (per-GFF genes.sqlite)."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from pathlib import Path

from helpers.curie import CurieError, format_curie, parse_curie
from helpers.cursor import (
    QueryError,
    decode_ann_genes_cursor,
    decode_gene_xrefs_cursor,
    encode_ann_genes_cursor,
    encode_gene_xrefs_cursor,
    filter_fingerprint,
    resolve_page_token,
)
from helpers.paths import per_gff_root, shard_genes_path
from helpers.pagination import DEFAULT_LIMIT, MAX_CURIES_GET, clamp_limit
from helpers.session import GeneCorpusUnavailable
from helpers.strand import format_strand
from services.hits import MATCH_ANY, parse_curie_list, parse_match

# helpers.curie bootstraps gene_corpus on sys.path.
from gene_corpus.namespaces.tiers import TIER_B  # noqa: E402

_SHARD_CACHE_SIZE = -16384  # 16 MiB
MAX_Q_CHARS = 64

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


def _gene_card(row: sqlite3.Row, *, annotation_id: str) -> dict:
    return {
        "annotation_id": annotation_id,
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
    *,
    annotation_id: str,
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
    return _gene_card(row, annotation_id=annotation_id)


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
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
) -> dict | None:
    """
    Paginated xref chips for one gene.

    Returns None when the gene is missing.
    """
    if not gene_exists(conn, local_id):
        return None

    limit = clamp_limit(limit)
    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_ns: str | None = None
    bound_acc: str | None = None
    if page_token is not None:
        bound_ns, bound_acc = decode_gene_xrefs_cursor(page_token)

    params: list = [int(local_id)]
    if direction == "next" and bound_ns is not None and bound_acc is not None:
        seek = (
            "AND (namespace > ? OR (namespace = ? AND accession > ?))"
        )
        params.extend([bound_ns, bound_ns, bound_acc])
        order = "ORDER BY namespace ASC, accession ASC"
    elif direction == "prev" and bound_ns is not None and bound_acc is not None:
        seek = (
            "AND (namespace < ? OR (namespace = ? AND accession < ?))"
        )
        params.extend([bound_ns, bound_ns, bound_acc])
        order = "ORDER BY namespace DESC, accession DESC"
    else:
        seek = ""
        order = "ORDER BY namespace ASC, accession ASC"

    rows = conn.execute(
        f"""
        SELECT namespace, accession, display, origin, via_level
        FROM gene_xref
        WHERE local_id = ?
        {seek}
        {order}
        LIMIT ?
        """,
        [*params, limit + 1],
    ).fetchall()

    going_prev = direction == "prev"
    has_more = len(rows) > limit
    rows = rows[:limit]
    if going_prev:
        rows = list(reversed(rows))

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

    next_token = None
    prev_token = None
    if results:
        first = results[0]
        last = results[-1]
        if going_prev:
            if has_more:
                prev_token = encode_gene_xrefs_cursor(
                    namespace=first["namespace"],
                    accession=first["accession"],
                )
            next_token = encode_gene_xrefs_cursor(
                namespace=last["namespace"],
                accession=last["accession"],
            )
        else:
            if has_more:
                next_token = encode_gene_xrefs_cursor(
                    namespace=last["namespace"],
                    accession=last["accession"],
                )
            if direction == "next":
                prev_token = encode_gene_xrefs_cursor(
                    namespace=first["namespace"],
                    accession=first["accession"],
                )

    return {
        "annotation_id": annotation_id,
        "local_id": int(local_id),
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }


def _escape_like(text: str) -> str:
    """Escape ``\\``, ``%``, and ``_`` for SQLite LIKE with ESCAPE '\\'."""
    return (
        text.replace("\\", "\\\\")
        .replace("%", "\\%")
        .replace("_", "\\_")
    )


def list_annotation_namespaces(
    conn: sqlite3.Connection,
    annotation_id: str,
) -> dict:
    """
    Namespaces present on one annotation shard (Tier A + Tier B).

    ``accession_count`` is the number of distinct accessions per namespace.
    Shape matches ``GET /namespaces``, plus ``annotation_id``.
    """
    by_ns: dict[str, int] = {}
    for row in conn.execute(
        """
        SELECT namespace, COUNT(*) AS accession_count
        FROM tier_a_counts
        GROUP BY namespace
        """
    ):
        by_ns[str(row["namespace"])] = int(row["accession_count"])

    tier_b = sorted(TIER_B)
    if tier_b:
        placeholders = ", ".join("?" for _ in tier_b)
        for row in conn.execute(
            f"""
            SELECT namespace, COUNT(DISTINCT accession) AS accession_count
            FROM gene_xref
            WHERE namespace IN ({placeholders})
            GROUP BY namespace
            """,
            tier_b,
        ):
            by_ns[str(row["namespace"])] = int(row["accession_count"])

    results = [
        {"namespace": ns, "accession_count": by_ns[ns]}
        for ns in sorted(by_ns)
    ]
    return {
        "annotation_id": annotation_id,
        "total": len(results),
        "limit": len(results),
        "results": results,
        "next": None,
        "previous": None,
    }


def _resolve_name_prefix(q: str | None) -> str | None:
    """Optional casefold prefix for ``primary_name`` (not a CURIE)."""
    if q is None:
        return None
    text = q.strip()
    if not text:
        return None
    if len(text) > MAX_Q_CHARS:
        raise QueryError(
            f"q accepts at most {MAX_Q_CHARS} characters",
            code="q_too_long",
            max=MAX_Q_CHARS,
            received=len(text),
        )
    return text.casefold()


def _parse_gene_curies(
    curies: str | Sequence[str] | None,
    *,
    match: str | None,
    max_curies: int,
) -> tuple[str, list[dict], dict[tuple[str, str], list[str]], list[str]]:
    """
    Parse optional ``curies`` + ``match`` for the genes list.

    Returns ``(match_mode, errors, by_key, display_ids)``.
    When ``curies`` is omitted, returns default match, empty errors/keys.
    """
    match_mode = parse_match(match)
    if curies is None:
        return match_mode, [], {}, []

    raw_list = parse_curie_list(curies, max_curies=max_curies)
    errors: list[dict] = []
    parsed: list[tuple[str, str]] = []
    for raw in raw_list:
        try:
            namespace, accession = parse_curie(raw)
            parsed.append((namespace, accession))
        except CurieError as exc:
            errors.append(
                {"curie": raw, "code": exc.code, "message": str(exc)}
            )

    by_key: dict[tuple[str, str], list[str]] = {}
    for namespace, accession in parsed:
        by_key.setdefault((namespace, accession), []).append(
            format_curie(namespace, accession)
        )
    display_ids = sorted(
        {format_curie(ns, acc) for ns, acc in by_key}
    )
    return match_mode, errors, by_key, display_ids


def list_annotation_genes(
    conn: sqlite3.Connection,
    annotation_id: str,
    *,
    q: str | None = None,
    curies: str | Sequence[str] | None = None,
    match: str | None = MATCH_ANY,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Paginated gene cards for one annotation shard.

    Optional ``q``: casefold prefix match on ``primary_name``.
    Optional ``curies``: xref filter (``match=any`` union / ``match=all``
    intersection). Both filters AND together when set.
    """
    limit = clamp_limit(limit)
    prefix = _resolve_name_prefix(q)
    match_mode, errors, by_key, display_ids = _parse_gene_curies(
        curies, match=match, max_curies=max_curies
    )
    curies_requested = curies is not None
    filter_f = filter_fingerprint(
        {
            "q": f"name:{prefix}" if prefix is not None else None,
            "curies": display_ids if curies_requested else None,
            "match": match_mode if curies_requested else None,
        }
    )

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_lid: int | None = None
    if page_token is not None:
        bound_lid = decode_ann_genes_cursor(page_token, filter_f=filter_f)

    # match=all with any bad prefix → empty page (same as /hits).
    if curies_requested and match_mode == "all" and errors:
        return {
            "annotation_id": annotation_id,
            "limit": limit,
            "results": [],
            "next": None,
            "previous": None,
            "errors": errors,
        }

    # match=any with only bad prefixes (no valid keys) → empty + errors.
    if curies_requested and not by_key:
        out = {
            "annotation_id": annotation_id,
            "limit": limit,
            "results": [],
            "next": None,
            "previous": None,
        }
        if errors:
            out["errors"] = errors
        return out

    params: list = []
    where_parts: list[str] = ["1=1"]

    if prefix is not None:
        params.append(_escape_like(prefix) + "%")
        where_parts.append(
            "g.primary_name IS NOT NULL "
            "AND LOWER(g.primary_name) LIKE ? ESCAPE '\\'"
        )

    if by_key:
        keys = list(by_key)
        in_tuples = ", ".join("(?, ?)" for _ in keys)
        for ns, acc in keys:
            params.extend([ns, acc])
        if match_mode == "all" and len(keys) > 1:
            params.append(len(keys))
            xref_sql = f"""
                SELECT x.local_id FROM gene_xref x
                WHERE (x.namespace, x.accession) IN ({in_tuples})
                GROUP BY x.local_id
                HAVING COUNT(*) = ?
            """
        else:
            xref_sql = f"""
                SELECT x.local_id FROM gene_xref x
                WHERE (x.namespace, x.accession) IN ({in_tuples})
            """
        where_parts.append(f"g.local_id IN ({xref_sql})")

    where_sql = " AND ".join(where_parts)

    if direction == "next" and bound_lid is not None:
        seek = "AND g.local_id > ?"
        params.append(bound_lid)
        order = "ORDER BY g.local_id ASC"
    elif direction == "prev" and bound_lid is not None:
        seek = "AND g.local_id < ?"
        params.append(bound_lid)
        order = "ORDER BY g.local_id DESC"
    else:
        seek = ""
        order = "ORDER BY g.local_id ASC"

    rows = conn.execute(
        f"""
        SELECT g.local_id, g.source_gene_id, g.feature_type, g.seqid,
               g.start, g.end, g.strand, g.biotype, g.primary_name,
               g.prose, g.prose_kind
        FROM gene g
        WHERE {where_sql}
        {seek}
        {order}
        LIMIT ?
        """,
        [*params, limit + 1],
    ).fetchall()

    going_prev = direction == "prev"
    has_more = len(rows) > limit
    rows = rows[:limit]
    if going_prev:
        rows = list(reversed(rows))

    results = [_gene_card(row, annotation_id=annotation_id) for row in rows]

    next_token = None
    prev_token = None
    if results:
        first_lid = int(results[0]["local_id"])
        last_lid = int(results[-1]["local_id"])
        if going_prev:
            if has_more:
                prev_token = encode_ann_genes_cursor(
                    local_id=first_lid, filter_f=filter_f
                )
            next_token = encode_ann_genes_cursor(
                local_id=last_lid, filter_f=filter_f
            )
        else:
            if has_more:
                next_token = encode_ann_genes_cursor(
                    local_id=last_lid, filter_f=filter_f
                )
            if direction == "next":
                prev_token = encode_ann_genes_cursor(
                    local_id=first_lid, filter_f=filter_f
                )

    out = {
        "annotation_id": annotation_id,
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }
    if errors:
        out["errors"] = errors
    return out
