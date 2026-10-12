"""Annotation presence index (xref_file) read paths."""
from __future__ import annotations

import sqlite3
from collections.abc import Sequence

from helpers.pagination import MAX_CURIES_GET
from helpers.source_codes import SOURCES_LEGEND, source_code
from services.hits import (
    MATCH_ALL,
    MATCH_ANY,
    _parse_hits_request,
    taxid_has_annotations,
)

ANNOTATIONS_INDEX_FIELDS = [
    "annotation_id",
    "annotation_key",
    "organism_name",
    "assembly_accession",
    "taxid",
    "source",
]

HITS_ANNOTATIONS_INDEX_FIELDS = [
    "annotation_id",
    "annotation_key",
    "taxid",
    "source",
]


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1",
        (name,),
    ).fetchone()
    return row is not None


def list_annotations_index(conn: sqlite3.Connection) -> dict:
    """
    Full placement catalogue: annotations with xref_file_status.placed = 1.

    Ordered by annotation_key (the client slot).
    """
    if not _table_exists(conn, "xref_file_status"):
        return {
            "n_without_tier_a": 0,
            "sources": dict(SOURCES_LEGEND),
            "fields": list(ANNOTATIONS_INDEX_FIELDS),
            "rows": [],
        }

    n_without = int(
        conn.execute(
            "SELECT COUNT(*) FROM xref_file_status WHERE placed = 0"
        ).fetchone()[0]
    )
    rows = conn.execute(
        """
        SELECT a.annotation_id, a.annotation_key, a.organism_name,
               a.assembly_accession, a.taxid, a.source_database
        FROM annotation a
        JOIN xref_file_status s ON s.annotation_key = a.annotation_key
        WHERE s.placed = 1
        ORDER BY a.annotation_key
        """
    ).fetchall()
    out_rows = [
        [
            row["annotation_id"],
            int(row["annotation_key"]),
            row["organism_name"],
            row["assembly_accession"],
            int(row["taxid"]),
            source_code(row["source_database"]),
        ]
        for row in rows
    ]
    return {
        "n_without_tier_a": n_without,
        "sources": dict(SOURCES_LEGEND),
        "fields": list(ANNOTATIONS_INDEX_FIELDS),
        "rows": out_rows,
    }


def _keys_any(
    conn: sqlite3.Connection,
    by_key: dict[tuple[str, str], list[str]],
    taxid_i: int | None,
) -> list[int]:
    """Annotation keys that have at least one of the CURIEs."""
    if not by_key:
        return []
    preds = " OR ".join(
        "(x.namespace = ? AND x.accession = ?)" for _ in by_key
    )
    params: list = []
    for namespace, accession in by_key:
        params.extend([namespace, accession])
    lineage = ""
    if taxid_i is not None:
        lineage = (
            "AND EXISTS ("
            "SELECT 1 FROM annotation_lineage l "
            "WHERE l.annotation_key = x.annotation_key AND l.taxid = ?"
            ")"
        )
        params.append(taxid_i)
    rows = conn.execute(
        f"""
        SELECT DISTINCT x.annotation_key AS annotation_key
        FROM xref_file x
        WHERE ({preds})
        {lineage}
        ORDER BY x.annotation_key
        """,
        params,
    ).fetchall()
    return [int(r["annotation_key"]) for r in rows]


def _keys_all(
    conn: sqlite3.Connection,
    by_key: dict[tuple[str, str], list[str]],
    taxid_i: int | None,
) -> list[int]:
    """Annotation keys that have every CURIE (annotation-level presence)."""
    if not by_key:
        return []
    keys_list = list(by_key.keys())
    # Start from the first posting list, then intersect.
    first_ns, first_acc = keys_list[0]
    params: list = [first_ns, first_acc]
    lineage = ""
    if taxid_i is not None:
        lineage = (
            "AND EXISTS ("
            "SELECT 1 FROM annotation_lineage l "
            "WHERE l.annotation_key = x.annotation_key AND l.taxid = ?"
            ")"
        )
        params.append(taxid_i)
    sql = f"""
        SELECT x.annotation_key AS annotation_key
        FROM xref_file x
        WHERE x.namespace = ? AND x.accession = ?
        {lineage}
    """
    for namespace, accession in keys_list[1:]:
        sql += """
        INTERSECT
        SELECT x.annotation_key AS annotation_key
        FROM xref_file x
        WHERE x.namespace = ? AND x.accession = ?
        """
        params.extend([namespace, accession])
        if taxid_i is not None:
            sql += """
            AND EXISTS (
                SELECT 1 FROM annotation_lineage l
                WHERE l.annotation_key = x.annotation_key AND l.taxid = ?
            )
            """
            params.append(taxid_i)
    sql += " ORDER BY 1"
    rows = conn.execute(sql, params).fetchall()
    return [int(r["annotation_key"]) for r in rows]


def _accession_present(
    conn: sqlite3.Connection, namespace: str, accession: str
) -> bool:
    row = conn.execute(
        """
        SELECT 1 FROM xref_file
        WHERE namespace = ? AND accession = ?
        LIMIT 1
        """,
        (namespace, accession),
    ).fetchone()
    return row is not None


def list_hit_annotations_index(
    conn: sqlite3.Connection,
    curies: str | Sequence[str],
    *,
    taxid: int | None = None,
    match: str | None = MATCH_ANY,
    max_curies: int = MAX_CURIES_GET,
) -> dict:
    """
    Full match set for CURIEs from xref_file (presence), ordered by key.

    Same match / taxid semantics as ``list_hit_annotations``, without paging
    and without n_genes / matched / organism fields.
    """
    match_mode, errors, by_key, _display_ids, _filter_f, taxid_i = (
        _parse_hits_request(
            curies, taxid=taxid, match=match, max_curies=max_curies
        )
    )

    def _envelope(*, rows: list) -> dict:
        # Same errors / taxid rules as /hits and /hits/annotations:
        # omit empty errors; echo taxid only when the request set it.
        out: dict = {
            "sources": dict(SOURCES_LEGEND),
            "fields": list(HITS_ANNOTATIONS_INDEX_FIELDS),
            "rows": rows,
        }
        if errors:
            out["errors"] = errors
        if taxid_i is not None:
            out["taxid"] = taxid_i
        return out

    empty = _envelope(rows=[])
    if not _table_exists(conn, "xref_file"):
        return empty

    if taxid_i is not None and not taxid_has_annotations(conn, taxid_i):
        return empty

    if match_mode == MATCH_ALL and errors:
        return empty

    if match_mode == MATCH_ALL and by_key:
        for namespace, accession in by_key:
            if not _accession_present(conn, namespace, accession):
                return empty

    if not by_key:
        return empty

    if match_mode == MATCH_ALL and len(by_key) > 1:
        order_keys = _keys_all(conn, by_key, taxid_i)
    else:
        order_keys = _keys_any(conn, by_key, taxid_i)

    if not order_keys:
        return empty

    placeholders = ",".join("?" for _ in order_keys)
    ann_rows = conn.execute(
        f"""
        SELECT annotation_key, annotation_id, taxid, source_database
        FROM annotation
        WHERE annotation_key IN ({placeholders})
        ORDER BY annotation_key
        """,
        order_keys,
    ).fetchall()
    out_rows = [
        [
            row["annotation_id"],
            int(row["annotation_key"]),
            int(row["taxid"]),
            source_code(row["source_database"]),
        ]
        for row in ann_rows
    ]
    return _envelope(rows=out_rows)
