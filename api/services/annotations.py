"""Annotation shard gene / xref queries (per-GFF genes.sqlite)."""
from __future__ import annotations

import sqlite3
from pathlib import Path

from helpers.curie import CurieError, parse_curie
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
from helpers.pagination import DEFAULT_LIMIT, MAX_LIMIT, clamp_limit
from helpers.session import GeneCorpusUnavailable
from helpers.strand import format_strand

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


def _resolve_gene_q(
    q: str | None,
) -> tuple[str, str, str] | tuple[str, str] | None:
    """
    Resolve optional ``q`` into a filter.

    Returns ``None`` (no filter), ``("curie", namespace, accession)``,
    or ``("name", prefix)``.
    """
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
    if ":" in text:
        try:
            namespace, accession = parse_curie(text)
        except CurieError as exc:
            raise QueryError(str(exc), code=exc.code) from exc
        return ("curie", namespace, accession)
    return ("name", text.casefold())


def list_annotation_genes(
    conn: sqlite3.Connection,
    annotation_id: str,
    *,
    q: str | None = None,
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
) -> dict:
    """
    Paginated gene cards for one annotation shard.

    Optional ``q``: CURIE (contains ``:``) filters via ``gene_xref``;
    otherwise casefold prefix match on ``primary_name``.
    """
    limit = clamp_limit(limit)
    resolved = _resolve_gene_q(q)
    if resolved is None:
        q_fp: str | None = None
    elif resolved[0] == "curie":
        q_fp = f"curie:{resolved[1]}:{resolved[2]}"
    else:
        q_fp = f"name:{resolved[1]}"
    filter_f = filter_fingerprint({"q": q_fp})

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_lid: int | None = None
    if page_token is not None:
        bound_lid = decode_ann_genes_cursor(page_token, filter_f=filter_f)

    params: list = []
    if resolved is None:
        base_from = "FROM gene g WHERE 1=1"
    elif resolved[0] == "curie":
        _kind, namespace, accession = resolved  # type: ignore[misc]
        params.extend([namespace, accession])
        base_from = """
            FROM gene g
            WHERE g.local_id IN (
                SELECT x.local_id FROM gene_xref x
                WHERE x.namespace = ? AND x.accession = ?
            )
        """
    else:
        _kind, prefix = resolved  # type: ignore[misc]
        params.append(_escape_like(prefix) + "%")
        base_from = """
            FROM gene g
            WHERE g.primary_name IS NOT NULL
              AND LOWER(g.primary_name) LIKE ? ESCAPE '\\'
        """

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
        {base_from}
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

    return {
        "annotation_id": annotation_id,
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }
