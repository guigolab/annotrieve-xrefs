#!/usr/bin/env python3
"""Audit species (taxid) coverage for ncbiortholog and eggnog in gene_corpus.

Uses the annotation presence index ``xref_file`` ⋈ ``annotation`` (same table
the finder uses for CURIE → annotations).

Reports per namespace:

* distinct taxids / annotations / accessions / xref_file rows
* breakdown by ``source_database``
* optional per-taxid table (counts + organism name)
* taxid set overlap between the two namespaces

Example::

    python3 -m gene_corpus.audit_ortho_taxid_coverage \\
        --db /data/annotrieve/gene_corpus/gene_corpus.sqlite

    python3 -m gene_corpus.audit_ortho_taxid_coverage \\
        --db /path/to/gene_corpus.sqlite --list-taxids --top 30
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
import time

DEFAULT_NAMESPACES = ("ncbiortholog", "eggnog")


def _pct(n: int, total: int) -> str:
    if total <= 0:
        return "n/a"
    return f"{100.0 * n / total:.2f}%"


def _corpus_baselines(conn: sqlite3.Connection) -> dict[str, int]:
    return {
        "annotations": int(
            conn.execute("SELECT COUNT(*) FROM annotation").fetchone()[0]
        ),
        "taxids": int(
            conn.execute(
                "SELECT COUNT(DISTINCT taxid) FROM annotation"
            ).fetchone()[0]
        ),
    }


def _namespace_summary(
    conn: sqlite3.Connection, namespace: str
) -> dict[str, int]:
    row = conn.execute(
        """
        SELECT
            COUNT(*) AS placed_rows,
            COUNT(DISTINCT x.accession) AS accessions,
            COUNT(DISTINCT x.annotation_key) AS annotations,
            COUNT(DISTINCT a.taxid) AS taxids
        FROM xref_file x
        JOIN annotation a ON a.annotation_key = x.annotation_key
        WHERE x.namespace = ?
        """,
        (namespace,),
    ).fetchone()
    meta_acc = int(
        conn.execute(
            "SELECT COUNT(*) FROM xref_meta WHERE namespace = ?",
            (namespace,),
        ).fetchone()[0]
    )
    return {
        "placed_rows": int(row[0]),
        "accessions_in_xref_file": int(row[1]),
        "annotations": int(row[2]),
        "taxids": int(row[3]),
        "accessions_in_meta": meta_acc,
    }


def _by_source(
    conn: sqlite3.Connection, namespace: str
) -> list[tuple[str, int, int]]:
    rows = conn.execute(
        """
        SELECT
            COALESCE(a.source_database, '') AS src,
            COUNT(DISTINCT a.taxid) AS taxids,
            COUNT(DISTINCT x.annotation_key) AS annotations
        FROM xref_file x
        JOIN annotation a ON a.annotation_key = x.annotation_key
        WHERE x.namespace = ?
        GROUP BY COALESCE(a.source_database, '')
        ORDER BY taxids DESC, annotations DESC, src
        """,
        (namespace,),
    ).fetchall()
    return [(str(r[0]), int(r[1]), int(r[2])) for r in rows]


def _taxid_rows(
    conn: sqlite3.Connection, namespace: str
) -> list[tuple[int, str, int, int, int]]:
    """(taxid, organism_name, annotations, accessions, placed_rows)."""
    rows = conn.execute(
        """
        SELECT
            a.taxid,
            COALESCE(
                MAX(NULLIF(a.organism_name, '')),
                ''
            ) AS organism,
            COUNT(DISTINCT x.annotation_key) AS annotations,
            COUNT(DISTINCT x.accession) AS accessions,
            COUNT(*) AS placed_rows
        FROM xref_file x
        JOIN annotation a ON a.annotation_key = x.annotation_key
        WHERE x.namespace = ?
        GROUP BY a.taxid
        ORDER BY annotations DESC, accessions DESC, a.taxid
        """,
        (namespace,),
    ).fetchall()
    return [
        (int(r[0]), str(r[1]), int(r[2]), int(r[3]), int(r[4]))
        for r in rows
    ]


def _taxid_set(conn: sqlite3.Connection, namespace: str) -> set[int]:
    rows = conn.execute(
        """
        SELECT DISTINCT a.taxid
        FROM xref_file x
        JOIN annotation a ON a.annotation_key = x.annotation_key
        WHERE x.namespace = ?
        """,
        (namespace,),
    ).fetchall()
    return {int(r[0]) for r in rows}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--db",
        required=True,
        help="path to gene_corpus.sqlite",
    )
    parser.add_argument(
        "--namespaces",
        nargs="+",
        default=list(DEFAULT_NAMESPACES),
        help=f"namespaces to audit (default: {' '.join(DEFAULT_NAMESPACES)})",
    )
    parser.add_argument(
        "--list-taxids",
        action="store_true",
        help="print per-taxid TSV for each namespace",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=0,
        help="with --list-taxids, only print top N taxids by annotations "
        "(0 = all)",
    )
    args = parser.parse_args(argv)
    namespaces = [ns.strip() for ns in args.namespaces if ns.strip()]
    if not namespaces:
        print("error: no namespaces", file=sys.stderr)
        return 2

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    try:
        t0 = time.monotonic()
        base = _corpus_baselines(conn)
        print("=== corpus baseline ===")
        print(f"annotations:  {base['annotations']:,}")
        print(f"distinct taxids: {base['taxids']:,}")
        print()

        taxid_sets: dict[str, set[int]] = {}

        for ns in namespaces:
            t1 = time.monotonic()
            print(f"=== {ns} ===")
            summary = _namespace_summary(conn, ns)
            print(f"xref_meta accessions:     {summary['accessions_in_meta']:,}")
            print(f"xref_file rows:           {summary['placed_rows']:,}")
            print(
                f"distinct accessions:      "
                f"{summary['accessions_in_xref_file']:,}"
            )
            print(f"annotations placed:       {summary['annotations']:,}")
            print(
                f"distinct taxids:          {summary['taxids']:,}  "
                f"({_pct(summary['taxids'], base['taxids'])} of corpus)"
            )
            print(
                f"annotations coverage:     "
                f"{_pct(summary['annotations'], base['annotations'])} "
                f"of corpus annotations"
            )

            print("by source_database (taxids / annotations):")
            for src, taxids, anns in _by_source(conn, ns):
                label = src or "(empty)"
                print(f"  {label:20s}  taxids={taxids:,}  annotations={anns:,}")

            taxid_sets[ns] = _taxid_set(conn, ns)
            print(f"query time: {time.monotonic() - t1:.1f}s")
            print()

            if args.list_taxids:
                rows = _taxid_rows(conn, ns)
                if args.top > 0:
                    rows = rows[: args.top]
                print(
                    f"# {ns} taxid\torganism\tannotations\t"
                    f"accessions\txref_file_rows"
                )
                for taxid, organism, anns, accs, placed in rows:
                    org = organism.replace("\t", " ").strip() or "-"
                    print(f"{taxid}\t{org}\t{anns}\t{accs}\t{placed}")
                print()

        if len(namespaces) >= 2:
            print("=== taxid set overlap ===")
            a, b = namespaces[0], namespaces[1]
            sa, sb = taxid_sets[a], taxid_sets[b]
            print(f"{a} taxids:           {len(sa):,}")
            print(f"{b} taxids:           {len(sb):,}")
            print(f"intersection:         {len(sa & sb):,}")
            print(f"only {a}:             {len(sa - sb):,}")
            print(f"only {b}:             {len(sb - sa):,}")
            print(
                f"union:                {len(sa | sb):,}  "
                f"({_pct(len(sa | sb), base['taxids'])} of corpus)"
            )
            if len(namespaces) > 2:
                inter_all = set.intersection(
                    *(taxid_sets[ns] for ns in namespaces)
                )
                print(
                    f"intersection all ({', '.join(namespaces)}): "
                    f"{len(inter_all):,}"
                )
            print()

        print(f"# total elapsed {time.monotonic() - t0:.1f}s", file=sys.stderr)
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
