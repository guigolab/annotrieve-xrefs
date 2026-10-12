"""Namespace and accession queries over gene_corpus.sqlite."""
from __future__ import annotations

import sqlite3

from helpers.curie import normalize_accession
from helpers.cursor import (
    QueryError,
    decode_accessions_cursor,
    decode_labels_cursor,
    encode_accessions_cursor,
    encode_labels_cursor,
    filter_fingerprint,
    resolve_page_token,
)
from helpers.namespace_names import namespace_name
from helpers.pagination import (
    DEFAULT_LIMIT,
    MAX_ACCESSIONS_FILTER,
    clamp_limit,
)

ACCESSION_SORTS = frozenset({"accession", "n_annotations", "n_genes"})
LABEL_SORTS = frozenset({"label", "accession", "n_annotations", "n_genes"})
MIN_LABEL_Q_LEN = 3


def _fts_match_query(q: str) -> str:
    """Quote *q* for FTS5 MATCH (trigram substring)."""
    return '"' + q.replace('"', '""') + '"'


def _namespaces_with_labels(conn: sqlite3.Connection) -> set[str]:
    """Namespaces that have at least one xref_label row."""
    try:
        rows = conn.execute(
            "SELECT DISTINCT namespace FROM xref_label"
        ).fetchall()
    except sqlite3.OperationalError:
        # Older/minimal DBs may lack the sidecar tables.
        return set()
    return {str(row["namespace"]) for row in rows}


def _namespace_row(
    namespace: str, accession_count: int, *, labels: bool
) -> dict:
    return {
        "namespace": namespace,
        "name": namespace_name(namespace),
        "accession_count": int(accession_count),
        "labels": bool(labels),
    }


def list_namespaces(
    conn: sqlite3.Connection,
    *,
    has_labels: bool | None = None,
) -> dict:
    labeled = _namespaces_with_labels(conn)
    rows = conn.execute(
        """
        SELECT namespace, accession_count
        FROM namespace_stats
        ORDER BY namespace
        """
    ).fetchall()
    results = [
        _namespace_row(
            row["namespace"],
            row["accession_count"],
            labels=row["namespace"] in labeled,
        )
        for row in rows
    ]
    if has_labels is not None:
        results = [row for row in results if row["labels"] is has_labels]
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
    labeled = _namespaces_with_labels(conn)
    return _namespace_row(
        row["namespace"],
        row["accession_count"],
        labels=row["namespace"] in labeled,
    )


def _accessions_page_sql(
    *,
    sort: str,
    sort_order: str,
    direction: str | None,
) -> tuple[str, str]:
    """Return (where_extra_sql, order_sql) for accession keyset pages."""
    descending = sort_order == "desc"
    # Public sort n_genes maps to SQL column n_loci.
    count_col = "n_loci" if sort == "n_genes" else "n_annotations"

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

    # sort == n_genes or n_annotations; stable tie-break accession ASC.
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


def parse_accession_filter(
    raw: str,
    *,
    namespace: str,
    max_accessions: int = MAX_ACCESSIONS_FILTER,
) -> list[str]:
    """Parse CSV accession filter; normalize and dedupe (order preserved)."""
    pieces: list[str] = []
    for part in str(raw).split(","):
        text = part.strip()
        if text:
            pieces.append(text)
    if len(pieces) > max_accessions:
        raise QueryError(
            f"accessions accepts at most {max_accessions} values",
            code="too_many_accessions",
            max=max_accessions,
            received=len(pieces),
        )
    out: list[str] = []
    seen: set[str] = set()
    for piece in pieces:
        canon = normalize_accession(namespace, piece)
        if canon not in seen:
            seen.add(canon)
            out.append(canon)
    return out


def get_accession(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
) -> dict | None:
    """O(1) xref_meta lookup; None if namespace or accession is unknown."""
    if get_namespace(conn, namespace) is None:
        return None
    canon = normalize_accession(namespace, accession)
    row = conn.execute(
        """
        SELECT namespace, accession, n_annotations, n_loci
        FROM xref_meta
        WHERE namespace = ? AND accession = ?
        """,
        (namespace, canon),
    ).fetchone()
    if row is None:
        return None
    return {
        "namespace": row["namespace"],
        "accession": row["accession"],
        "n_annotations": int(row["n_annotations"]),
        "n_genes": int(row["n_loci"]),
    }


def list_accessions(
    conn: sqlite3.Connection,
    namespace: str,
    *,
    sort: str = "n_genes",
    sort_order: str = "desc",
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
    accessions: str | None = None,
) -> dict | None:
    stats = get_namespace(conn, namespace)
    if stats is None:
        return None

    column = sort if sort in ACCESSION_SORTS else "n_genes"
    order = "asc" if (sort_order or "").lower() == "asc" else "desc"
    limit = clamp_limit(limit)

    filter_vals: list[str] | None = None
    if accessions is not None and str(accessions).strip() != "":
        filter_vals = parse_accession_filter(accessions, namespace=namespace)

    filter_f = filter_fingerprint(
        {
            "sort": column,
            "sort_order": order,
            "accessions": list(filter_vals) if filter_vals is not None else [],
        }
    )

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_n: int | None = None
    bound_a: str | None = None
    if page_token is not None:
        bound_n, bound_a = decode_accessions_cursor(
            page_token,
            sort=column,
            sort_order=order,
            filter_f=filter_f,
        )

    where_extra, order_sql = _accessions_page_sql(
        sort=column, sort_order=order, direction=direction
    )
    params: list = [namespace]
    filter_sql = ""
    if filter_vals is not None:
        if not filter_vals:
            return {
                "namespace": namespace,
                "total": 0,
                "limit": limit,
                "results": [],
                "next": None,
                "previous": None,
            }
        placeholders = ",".join("?" for _ in filter_vals)
        filter_sql = f"AND accession IN ({placeholders})"
        params.extend(filter_vals)

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
        {filter_sql}
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
            "n_genes": int(row["n_loci"]),
        }
        for row in rows
    ]

    def _sort_value(item: dict) -> int:
        if column == "accession":
            return 0
        return int(item["n_genes" if column == "n_genes" else "n_annotations"])

    next_token = None
    prev_token = None
    if results:
        first = results[0]
        last = results[-1]
        if going_prev:
            if has_more_in_scan:
                prev_token = encode_accessions_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(first),
                    accession=first["accession"],
                    filter_f=filter_f,
                )
            next_token = encode_accessions_cursor(
                sort=column,
                sort_order=order,
                sort_value=_sort_value(last),
                accession=last["accession"],
                filter_f=filter_f,
            )
        else:
            if has_more_in_scan:
                next_token = encode_accessions_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(last),
                    accession=last["accession"],
                    filter_f=filter_f,
                )
            if direction == "next":
                prev_token = encode_accessions_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_sort_value(first),
                    accession=first["accession"],
                    filter_f=filter_f,
                )

    if filter_vals is not None:
        total_row = conn.execute(
            f"""
            SELECT COUNT(*) AS n
            FROM xref_meta
            WHERE namespace = ?
            AND accession IN ({",".join("?" for _ in filter_vals)})
            """,
            [namespace, *filter_vals],
        ).fetchone()
        total = int(total_row["n"])
    else:
        total = stats["accession_count"]

    return {
        "namespace": namespace,
        "total": total,
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }


def _labels_meta_join(namespace: str) -> tuple[str, list]:
    """
    Return ``(join_sql, params)`` after ``FROM xref_label l``.

    Pfam labels are keyed unversioned (``PF00001``) while ``xref_meta`` may
    keep ``PF00001.15``; aggregate version spellings under the unversioned key.
    """
    if namespace == "pfam":
        return (
            """
        LEFT JOIN (
            SELECT
                CASE
                    WHEN upper(accession) GLOB 'PF[0-9]*'
                         AND instr(accession, '.') > 0
                    THEN upper(
                        substr(accession, 1, instr(accession, '.') - 1)
                    )
                    ELSE upper(accession)
                END AS accession,
                SUM(n_annotations) AS n_annotations,
                SUM(n_loci) AS n_loci
            FROM xref_meta
            WHERE namespace = ?
            GROUP BY 1
        ) m ON m.accession = l.accession
            """,
            [namespace],
        )
    return (
        """
        LEFT JOIN xref_meta m
            ON m.namespace = l.namespace AND m.accession = l.accession
        """,
        [],
    )


def _labels_page_sql(
    *,
    sort: str,
    sort_order: str,
    direction: str | None,
) -> tuple[str, str]:
    """Return (where_extra_sql, order_sql) for labels keyset pages."""
    descending = sort_order == "desc"
    count_col = "n_loci" if sort == "n_genes" else "n_annotations"
    # Qualified columns from xref_label l LEFT JOIN xref_meta m.
    acc = "l.accession"
    lab = "l.label"
    # COALESCE so NULL meta sorts as 0.
    cnt = f"COALESCE(m.{count_col}, 0)"

    if sort == "accession":
        if direction is None:
            order = f"{acc} DESC" if descending else f"{acc} ASC"
            return "", order
        if direction == "next":
            if descending:
                return f"AND {acc} < ?", f"{acc} DESC"
            return f"AND {acc} > ?", f"{acc} ASC"
        if descending:
            return f"AND {acc} > ?", f"{acc} ASC"
        return f"AND {acc} < ?", f"{acc} DESC"

    if sort == "label":
        if direction is None:
            if descending:
                return "", f"{lab} DESC, {acc} ASC"
            return "", f"{lab} ASC, {acc} ASC"
        if direction == "next":
            if descending:
                return (
                    f"AND ({lab} < ? OR ({lab} = ? AND {acc} > ?))",
                    f"{lab} DESC, {acc} ASC",
                )
            return (
                f"AND ({lab} > ? OR ({lab} = ? AND {acc} > ?))",
                f"{lab} ASC, {acc} ASC",
            )
        if descending:
            return (
                f"AND ({lab} > ? OR ({lab} = ? AND {acc} < ?))",
                f"{lab} ASC, {acc} DESC",
            )
        return (
            f"AND ({lab} < ? OR ({lab} = ? AND {acc} < ?))",
            f"{lab} DESC, {acc} DESC",
        )

    # sort == n_genes or n_annotations; stable tie-break accession ASC.
    if direction is None:
        if descending:
            return "", f"{cnt} DESC, {acc} ASC"
        return "", f"{cnt} ASC, {acc} ASC"

    if direction == "next":
        if descending:
            return (
                f"AND ({cnt} < ? OR ({cnt} = ? AND {acc} > ?))",
                f"{cnt} DESC, {acc} ASC",
            )
        return (
            f"AND ({cnt} > ? OR ({cnt} = ? AND {acc} > ?))",
            f"{cnt} ASC, {acc} ASC",
        )

    if descending:
        return (
            f"AND ({cnt} > ? OR ({cnt} = ? AND {acc} < ?))",
            f"{cnt} ASC, {acc} DESC",
        )
    return (
        f"AND ({cnt} < ? OR ({cnt} = ? AND {acc} < ?))",
        f"{cnt} DESC, {acc} DESC",
    )


def list_labels(
    conn: sqlite3.Connection,
    namespace: str,
    *,
    q: str | None = None,
    sort: str = "label",
    sort_order: str = "asc",
    next: str | None = None,
    previous: str | None = None,
    limit: int = DEFAULT_LIMIT,
) -> dict:
    """
    Paginated labels for a namespace that already has ``labels: true``.

    Caller must verify the namespace exists and has label rows.
    """
    column = sort if sort in LABEL_SORTS else "label"
    order = "asc" if (sort_order or "").lower() == "asc" else "desc"
    limit = clamp_limit(limit)

    q_text: str | None = None
    if q is not None and str(q).strip() != "":
        q_text = str(q).strip()
        if len(q_text) < MIN_LABEL_Q_LEN:
            raise QueryError(
                f"q must be at least {MIN_LABEL_Q_LEN} characters",
                code="invalid_q",
                min_length=MIN_LABEL_Q_LEN,
            )

    filter_f = filter_fingerprint(
        {
            "sort": column,
            "sort_order": order,
            "q": q_text or "",
        }
    )

    direction, page_token = resolve_page_token(next=next, previous=previous)
    bound_key: int | str | None = None
    bound_a: str | None = None
    if page_token is not None:
        bound_key, bound_a = decode_labels_cursor(
            page_token,
            sort=column,
            sort_order=order,
            filter_f=filter_f,
        )

    where_extra, order_sql = _labels_page_sql(
        sort=column, sort_order=order, direction=direction
    )

    meta_join, meta_params = _labels_meta_join(namespace)
    fts_sql = ""
    params: list = [*meta_params, namespace]
    if q_text is not None:
        fts_sql = (
            "AND l.accession IN ("
            "SELECT DISTINCT accession FROM xref_label_fts "
            "WHERE namespace = ? AND text MATCH ?"
            ")"
        )
        params.extend([namespace, _fts_match_query(q_text)])

    if direction is not None:
        assert bound_a is not None and bound_key is not None
        if column == "accession":
            params.append(bound_a)
        else:
            # label / count sorts: (key, key, accession)
            params.extend([bound_key, bound_key, bound_a])

    rows = conn.execute(
        f"""
        SELECT
            l.accession AS accession,
            l.label AS label,
            COALESCE(m.n_annotations, 0) AS n_annotations,
            COALESCE(m.n_loci, 0) AS n_loci
        FROM xref_label l
        {meta_join}
        WHERE l.namespace = ?
        {fts_sql}
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
            "label": row["label"],
            "n_annotations": int(row["n_annotations"]),
            "n_genes": int(row["n_loci"]),
        }
        for row in rows
    ]

    def _cursor_sort_value(item: dict) -> int | str:
        if column == "label":
            return item["label"]
        if column == "accession":
            return 0
        return int(item["n_genes" if column == "n_genes" else "n_annotations"])

    next_token = None
    prev_token = None
    if results:
        first = results[0]
        last = results[-1]
        if going_prev:
            if has_more_in_scan:
                prev_token = encode_labels_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_cursor_sort_value(first),
                    accession=first["accession"],
                    filter_f=filter_f,
                )
            next_token = encode_labels_cursor(
                sort=column,
                sort_order=order,
                sort_value=_cursor_sort_value(last),
                accession=last["accession"],
                filter_f=filter_f,
            )
        else:
            if has_more_in_scan:
                next_token = encode_labels_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_cursor_sort_value(last),
                    accession=last["accession"],
                    filter_f=filter_f,
                )
            if direction == "next":
                prev_token = encode_labels_cursor(
                    sort=column,
                    sort_order=order,
                    sort_value=_cursor_sort_value(first),
                    accession=first["accession"],
                    filter_f=filter_f,
                )

    if q_text is not None:
        total_row = conn.execute(
            """
            SELECT COUNT(DISTINCT accession) AS n
            FROM xref_label_fts
            WHERE namespace = ? AND text MATCH ?
            """,
            (namespace, _fts_match_query(q_text)),
        ).fetchone()
        total = int(total_row["n"])
    else:
        total_row = conn.execute(
            "SELECT COUNT(*) AS n FROM xref_label WHERE namespace = ?",
            (namespace,),
        ).fetchone()
        total = int(total_row["n"])

    return {
        "namespace": namespace,
        "total": total,
        "limit": limit,
        "results": results,
        "next": next_token,
        "previous": prev_token,
    }


