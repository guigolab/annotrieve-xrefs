"""Namespace directory, /matches facet, and scoped /terms queries."""
from __future__ import annotations

import re
import sqlite3

from helpers.curie import (
    CurieError,
    format_curie,
    normalize_accession,
    parse_curie,
)
from helpers.cursor import QueryError
from helpers.namespace_names import namespace_name
from helpers.pagination import (
    DEFAULT_LIMIT,
    MIN_Q_LEN,
    clamp_limit,
    clamp_offset,
    page_envelope,
)

from gene_corpus.labels.providers.pfam import (  # noqa: E402
    normalize_pfam_accession,
)
from gene_corpus.namespaces.tiers import (  # noqa: E402
    LABELED_SEARCH_NAMESPACES,
    SEARCHABLE_NAMESPACES,
)

_DIGITS_RE = re.compile(r"^\d+$")


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


def _is_searchable(namespace: str) -> bool:
    return namespace in SEARCHABLE_NAMESPACES


def list_namespaces(
    conn: sqlite3.Connection,
    *,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
) -> dict:
    """Allowlisted namespaces present in namespace_stats (fixed shape)."""
    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    labeled = _namespaces_with_labels(conn)
    rows = conn.execute(
        """
        SELECT namespace, accession_count
        FROM namespace_stats
        ORDER BY namespace
        """
    ).fetchall()
    all_rows = [
        _namespace_row(
            row["namespace"],
            row["accession_count"],
            labels=row["namespace"] in labeled,
        )
        for row in rows
        if _is_searchable(str(row["namespace"]))
    ]
    total = len(all_rows)
    page = all_rows[offset : offset + limit]
    return page_envelope(
        total=total, limit=limit, offset=offset, results=page
    )


def get_namespace(conn: sqlite3.Connection, namespace: str) -> dict | None:
    if not _is_searchable(namespace):
        return None
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


def _meta_exact_count(
    conn: sqlite3.Connection, namespace: str, accession: str
) -> int:
    row = conn.execute(
        """
        SELECT 1 FROM xref_meta
        WHERE namespace = ? AND accession = ?
        LIMIT 1
        """,
        (namespace, accession),
    ).fetchone()
    return 1 if row is not None else 0


def _meta_prefix_count(
    conn: sqlite3.Connection, namespace: str, prefix: str
) -> int:
    row = conn.execute(
        """
        SELECT COUNT(*) AS n FROM xref_meta
        WHERE namespace = ? AND accession >= ? AND accession < ?
        """,
        (namespace, prefix, prefix + "\uffff"),
    ).fetchone()
    return int(row["n"])


def _labeled_accession_prefix(namespace: str, q: str) -> str:
    """Prefix for xref_label.accession matching (Pfam strips .version)."""
    prefix = normalize_accession(namespace, q)
    if namespace == "pfam":
        stripped = normalize_pfam_accession(prefix)
        if stripped is not None:
            return stripped
    return prefix


def _labeled_match_count(
    conn: sqlite3.Connection, namespace: str, q: str
) -> int:
    """Distinct label accessions matching FTS text or accession prefix."""
    prefix = _labeled_accession_prefix(namespace, q)
    try:
        row = conn.execute(
            """
            SELECT COUNT(*) AS n FROM (
                SELECT accession
                FROM xref_label_fts
                WHERE namespace = ? AND text MATCH ?
                UNION
                SELECT accession
                FROM xref_label
                WHERE namespace = ?
                  AND accession >= ?
                  AND accession < ?
            )
            """,
            (
                namespace,
                _fts_match_query(q),
                namespace,
                prefix,
                prefix + "\uffff",
            ),
        ).fetchone()
    except sqlite3.OperationalError:
        # FTS unavailable — accession prefix only
        row = conn.execute(
            """
            SELECT COUNT(*) AS n
            FROM xref_label
            WHERE namespace = ?
              AND accession >= ?
              AND accession < ?
            """,
            (namespace, prefix, prefix + "\uffff"),
        ).fetchone()
    return int(row["n"])


def _match_count_for_namespace(
    conn: sqlite3.Connection, namespace: str, q: str
) -> int:
    """Exact hit count for one searchable namespace and query string."""
    if namespace in LABELED_SEARCH_NAMESPACES:
        return _labeled_match_count(conn, namespace, q)
    if namespace in {"symbol", "alias"}:
        prefix = normalize_accession(namespace, q)
        return _meta_prefix_count(conn, namespace, prefix)
    if namespace in {"geneid", "ncbiortholog"}:
        canon = normalize_accession(namespace, q)
        return _meta_exact_count(conn, namespace, canon)
    # cdd, eggnog, …
    canon = normalize_accession(namespace, q)
    return _meta_prefix_count(conn, namespace, canon)


def list_namespace_matches(
    conn: sqlite3.Connection,
    q: str,
    *,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
) -> dict:
    """
    Facet: which searchable namespaces match *q* (always includes match_count).
    """
    text = (q or "").strip()
    if not text:
        raise QueryError("q is required", code="q_required")

    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    labeled = _namespaces_with_labels(conn)

    # CURIE → single namespace
    if ":" in text:
        try:
            namespace, accession = parse_curie(text)
        except CurieError as exc:
            raise QueryError(str(exc), code=exc.code) from exc
        count = _meta_exact_count(conn, namespace, accession)
        results = []
        if count > 0:
            stats = get_namespace(conn, namespace)
            if stats is not None:
                results.append({**stats, "match_count": count})
        return page_envelope(
            total=len(results),
            limit=limit,
            offset=offset,
            results=results[offset : offset + limit],
        )

    # Digits → geneid (and ncbiortholog if present as exact)
    if _DIGITS_RE.fullmatch(text):
        hits: list[dict] = []
        for ns in ("geneid", "ncbiortholog"):
            count = _meta_exact_count(conn, ns, text)
            if count > 0:
                stats = get_namespace(conn, ns)
                if stats is not None:
                    hits.append({**stats, "match_count": count})
        hits.sort(key=lambda r: r["namespace"])
        return page_envelope(
            total=len(hits),
            limit=limit,
            offset=offset,
            results=hits[offset : offset + limit],
        )

    if len(text) < MIN_Q_LEN:
        raise QueryError(
            f"q must be at least {MIN_Q_LEN} characters",
            code="q_too_short",
            min_length=MIN_Q_LEN,
        )

    # Directory rows for allowlisted ns present in stats
    stats_rows = {
        str(row["namespace"]): int(row["accession_count"])
        for row in conn.execute(
            "SELECT namespace, accession_count FROM namespace_stats"
        )
        if _is_searchable(str(row["namespace"]))
    }

    hits = []
    for namespace in sorted(stats_rows):
        count = _match_count_for_namespace(conn, namespace, text)
        if count <= 0:
            continue
        hits.append(
            {
                **_namespace_row(
                    namespace,
                    stats_rows[namespace],
                    labels=namespace in labeled,
                ),
                "match_count": count,
            }
        )

    return page_envelope(
        total=len(hits),
        limit=limit,
        offset=offset,
        results=hits[offset : offset + limit],
    )


def _labels_meta_join(namespace: str) -> tuple[str, list]:
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


def _term_row(
    namespace: str,
    accession: str,
    *,
    label: str | None,
    n_annotations: int,
    n_genes: int,
) -> dict:
    return {
        "curie": format_curie(namespace, accession),
        "namespace": namespace,
        "accession": accession,
        "label": label,
        "n_annotations": int(n_annotations),
        "n_genes": int(n_genes),
    }


def get_term(
    conn: sqlite3.Connection,
    namespace: str,
    accession: str,
) -> dict | None:
    """One term card from xref_meta (+ label when present)."""
    if not _is_searchable(namespace):
        return None
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

    label = None
    try:
        lab = conn.execute(
            """
            SELECT label FROM xref_label
            WHERE namespace = ? AND accession = ?
            """,
            (namespace, row["accession"]),
        ).fetchone()
        if lab is not None:
            label = str(lab["label"])
    except sqlite3.OperationalError:
        pass

    return _term_row(
        namespace,
        str(row["accession"]),
        label=label,
        n_annotations=int(row["n_annotations"]),
        n_genes=int(row["n_loci"]),
    )


def list_terms(
    conn: sqlite3.Connection,
    namespace: str,
    *,
    q: str | None = None,
    limit: int = DEFAULT_LIMIT,
    offset: int = 0,
) -> dict | None:
    """
    Paged terms in one searchable namespace.

    Returns None when the namespace is unknown / not searchable.
    """
    stats = get_namespace(conn, namespace)
    if stats is None:
        return None

    limit = clamp_limit(limit)
    offset = clamp_offset(offset)
    q_text = q.strip() if q is not None and str(q).strip() else None

    # Labeled namespaces: FTS or full label browse
    if namespace in LABELED_SEARCH_NAMESPACES and stats["labels"]:
        if q_text is not None and len(q_text) < MIN_Q_LEN:
            raise QueryError(
                f"q must be at least {MIN_Q_LEN} characters",
                code="q_too_short",
                min_length=MIN_Q_LEN,
            )
        meta_join, meta_params = _labels_meta_join(namespace)
        q_filter = ""
        params: list = [*meta_params, namespace]
        if q_text is not None:
            prefix = _labeled_accession_prefix(namespace, q_text)
            q_filter = """
                AND l.accession IN (
                    SELECT accession
                    FROM xref_label_fts
                    WHERE namespace = ? AND text MATCH ?
                    UNION
                    SELECT accession
                    FROM xref_label
                    WHERE namespace = ?
                      AND accession >= ?
                      AND accession < ?
                )
            """
            params.extend(
                [
                    namespace,
                    _fts_match_query(q_text),
                    namespace,
                    prefix,
                    prefix + "\uffff",
                ]
            )
            total = _labeled_match_count(conn, namespace, q_text)
        else:
            total_row = conn.execute(
                "SELECT COUNT(*) AS n FROM xref_label WHERE namespace = ?",
                (namespace,),
            ).fetchone()
            total = int(total_row["n"])

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
            {q_filter}
            ORDER BY l.label ASC, l.accession ASC
            LIMIT ? OFFSET ?
            """,
            [*params, limit, offset],
        ).fetchall()
        results = [
            _term_row(
                namespace,
                str(row["accession"]),
                label=str(row["label"]),
                n_annotations=int(row["n_annotations"]),
                n_genes=int(row["n_loci"]),
            )
            for row in rows
        ]
        return page_envelope(
            total=total,
            limit=limit,
            offset=offset,
            results=results,
            namespace=namespace,
        )

    # Identity / id-like: require q
    if q_text is None:
        raise QueryError(
            "q is required for this namespace",
            code="q_required",
        )
    if namespace in {"geneid", "ncbiortholog"}:
        # exact only
        term = get_term(conn, namespace, q_text)
        results = [term] if term is not None else []
        return page_envelope(
            total=len(results),
            limit=limit,
            offset=offset,
            results=results[offset : offset + limit],
            namespace=namespace,
        )

    if len(q_text) < MIN_Q_LEN and namespace in {"symbol", "alias"}:
        raise QueryError(
            f"q must be at least {MIN_Q_LEN} characters",
            code="q_too_short",
            min_length=MIN_Q_LEN,
        )

    prefix = normalize_accession(namespace, q_text)
    total = _meta_prefix_count(conn, namespace, prefix)
    rows = conn.execute(
        """
        SELECT accession, n_annotations, n_loci
        FROM xref_meta
        WHERE namespace = ? AND accession >= ? AND accession < ?
        ORDER BY n_annotations DESC, accession ASC
        LIMIT ? OFFSET ?
        """,
        (namespace, prefix, prefix + "\uffff", limit, offset),
    ).fetchall()
    results = [
        _term_row(
            namespace,
            str(row["accession"]),
            label=None,
            n_annotations=int(row["n_annotations"]),
            n_genes=int(row["n_loci"]),
        )
        for row in rows
    ]
    return page_envelope(
        total=total,
        limit=limit,
        offset=offset,
        results=results,
        namespace=namespace,
    )
