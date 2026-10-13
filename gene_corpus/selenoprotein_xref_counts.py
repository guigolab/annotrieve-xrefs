#!/usr/bin/env python3
"""Count xref_file presence for symbol/alias genes by annotation source_database.

Example::

    python3 -m gene_corpus.selenoprotein_xref_counts \\
        --db /data/annotrieve/gene_corpus/gene_corpus.sqlite
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import defaultdict

# Duplicates in the source list are collapsed (order preserved).
GENES = [
    "GPX1",
    "GPX2",
    "GPX3",
    "GPX4",
    "GPX6",
    "DIO1",
    "DIO2",
    "DIO3",
    "SELENON",
    "TXNRD1",
    "TXNRD2",
    "TXNRD3",
    "SELENOK",
    "SELENOP",
    "SELENOH",
    "SELENOW",
    "SELENOF",
    "SELENOM",
    "SELENOI",
    "SELENOS",
    "SELENOT",
    "SELENOV",
    "SELENOO",
    "SEPHS2",
    "MSRB1",
    "SELENOJ",
    "SELENOU",
    "SELENOL",
    "SELENOE",
]

SOURCES = ("Ensembl", "GenBank", "RefSeq")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        required=True,
        help="path to gene_corpus.sqlite",
    )
    args = parser.parse_args(argv)

    # Harvest stores symbol/alias casefolded.
    accessions = list(
        dict.fromkeys(g.strip().casefold() for g in GENES if g.strip())
    )
    placeholders = ",".join("?" for _ in accessions)

    sql = f"""
    SELECT
        x.namespace,
        x.accession,
        lower(a.source_database) AS src,
        COUNT(*) AS n
    FROM xref_file x
    JOIN annotation a ON a.annotation_key = x.annotation_key
    WHERE x.namespace IN ('symbol', 'alias')
      AND x.accession IN ({placeholders})
    GROUP BY x.namespace, x.accession, lower(a.source_database)
    ORDER BY x.accession, x.namespace, src
    """

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(sql, accessions).fetchall()
    finally:
        conn.close()

    # gene -> namespace -> source -> count
    tallies: dict[str, dict[str, dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )
    for r in rows:
        src = str(r["src"] or "")
        label = next((s for s in SOURCES if s.casefold() == src), None)
        if label is None:
            label = f"other:{src}" if src else "other"
        tallies[r["accession"]][r["namespace"]][label] = int(r["n"])

    cols = ["gene", "ns"] + list(SOURCES) + ["other", "total"]
    print("\t".join(cols))

    for gene in accessions:
        for ns in ("symbol", "alias"):
            by_src = tallies.get(gene, {}).get(ns, {})
            if not by_src:
                continue
            ens = by_src.get("Ensembl", 0)
            gb = by_src.get("GenBank", 0)
            rs = by_src.get("RefSeq", 0)
            other = sum(v for k, v in by_src.items() if k not in SOURCES)
            total = ens + gb + rs + other
            print(f"{gene}\t{ns}\t{ens}\t{gb}\t{rs}\t{other}\t{total}")

        if gene not in tallies:
            print(f"{gene}\t-\t0\t0\t0\t0\t0")

    print(
        "# note: symbol+alias rows are separate; "
        "an annotation may count in both",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
