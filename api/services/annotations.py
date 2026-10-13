"""Annotation shard gene / xref queries (per-GFF genes.sqlite)."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence
from pathlib import Path

from helpers.curie import CurieError, format_curie, parse_curie
from helpers.cursor import QueryError
from helpers.match import MATCH_ANY, parse_curie_list, parse_match
from helpers.paths import per_gff_root, shard_genes_path
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_CURIES_GET,
    clamp_limit,
    clamp_offset,
    page_envelope,
)
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
    return page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
        annotation_id=annotation_id,
        local_id=int(local_id),
    )


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
    All namespaces present on one annotation shard (from gene_xref).

    ``accession_count`` is the number of distinct accessions per namespace.
    """
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
        annotation_id=annotation_id,
    )


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
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Paginated gene cards for one annotation shard.

    Optional ``q``: casefold prefix match on ``primary_name``.
    Optional ``curies``: xref filter (``match=any`` / ``match=all``).
    """
    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    prefix = _resolve_name_prefix(q)
    match_mode, errors, by_key, _display_ids = _parse_gene_curies(
        curies, match=match, max_curies=max_curies
    )
    curies_requested = curies is not None

    empty = page_envelope(
        total=0,
        limit=limit,
        offset=offset,
        results=[],
        annotation_id=annotation_id,
    )
    if errors:
        empty["errors"] = errors

    if curies_requested and match_mode == "all" and errors:
        return empty

    if curies_requested and not by_key:
        return empty

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

    results = [_gene_card(row, annotation_id=annotation_id) for row in rows]
    out = page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
        annotation_id=annotation_id,
    )
    if errors:
        out["errors"] = errors
    return out
