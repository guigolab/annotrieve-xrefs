#!/usr/bin/env python3
"""Audit whether ncbiortholog accessions relate to geneid in gene_corpus.

Reports two relationships:

1. **Accession equality** (xref_meta): the ortholog group id also exists as a
   ``geneid`` accession (common when NCBI uses a GeneID as the group id).
2. **Same-locus co-occurrence** (gene_hit): at least one gene that carries
   ``ncbiortholog:X`` also carries some ``geneid`` (any accession) on the same
   ``(annotation_key, local_id)``.

Example::

    python3 -m gene_corpus.audit_ncbiortholog_geneid \\
        --db /data/annotrieve/gene_corpus/gene_corpus.sqlite

    python3 -m gene_corpus.audit_ncbiortholog_geneid \\
        --db /path/to/gene_corpus.sqlite --sample 20 --skip-loci
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
import time
from collections.abc import Iterable


def _pct(n: int, total: int) -> str:
    if total <= 0:
        return "n/a"
    return f"{100.0 * n / total:.2f}%"


def _print_sample(title: str, rows: Iterable[str], limit: int) -> None:
    items = list(rows)
    if not items:
        print(f"{title}: (none)")
        return
    print(f"{title} (showing {min(limit, len(items))} of {len(items)}):")
    for acc in items[:limit]:
        print(f"  {acc}")


def audit_meta(
    conn: sqlite3.Connection, *, sample: int
) -> tuple[int, int, list[str]]:
    """Return (total_ortho, with_equal_geneid, sample_missing)."""
    total = int(
        conn.execute(
            "SELECT COUNT(*) FROM xref_meta WHERE namespace = 'ncbiortholog'"
        ).fetchone()[0]
    )
    with_geneid = int(
        conn.execute(
            """
            SELECT COUNT(*)
            FROM xref_meta o
            WHERE o.namespace = 'ncbiortholog'
              AND EXISTS (
                  SELECT 1
                  FROM xref_meta g
                  WHERE g.namespace = 'geneid'
                    AND g.accession = o.accession
              )
            """
        ).fetchone()[0]
    )
    missing_rows = conn.execute(
        """
        SELECT o.accession
        FROM xref_meta o
        WHERE o.namespace = 'ncbiortholog'
          AND NOT EXISTS (
              SELECT 1
              FROM xref_meta g
              WHERE g.namespace = 'geneid'
                AND g.accession = o.accession
          )
        ORDER BY o.accession
        LIMIT ?
        """,
        (sample,),
    ).fetchall()
    missing_sample = [str(r[0]) for r in missing_rows]
    return total, with_geneid, missing_sample


def audit_loci(
    conn: sqlite3.Connection, *, sample: int
) -> dict[str, int | list[str]]:
    """
    Stream gene_hit and measure same-locus co-occurrence.

    Memory scales with distinct ncbiortholog loci and accessions (~ortho size),
    not with full geneid placement count.
    """
    t0 = time.monotonic()
    # accession -> set of loci that carry it; also all ortho loci
    ortho_loci: dict[tuple[int, int], set[str]] = {}
    placement_rows = 0
    for annotation_key, local_id, accession in conn.execute(
        """
        SELECT annotation_key, local_id, accession
        FROM gene_hit
        WHERE namespace = 'ncbiortholog'
        """
    ):
        placement_rows += 1
        key = (int(annotation_key), int(local_id))
        ortho_loci.setdefault(key, set()).add(str(accession))
        if placement_rows % 500_000 == 0:
            print(
                f"# loci pass1: {placement_rows:,} ncbiortholog gene_hit rows…",
                file=sys.stderr,
                flush=True,
            )

    distinct_accessions = {
        acc for accs in ortho_loci.values() for acc in accs
    }
    print(
        f"# loci pass1 done: {placement_rows:,} rows, "
        f"{len(ortho_loci):,} loci, {len(distinct_accessions):,} accessions "
        f"({time.monotonic() - t0:.1f}s)",
        file=sys.stderr,
        flush=True,
    )

    loci_with_geneid: set[tuple[int, int]] = set()
    geneid_rows_scanned = 0
    t1 = time.monotonic()
    # Only need to mark ortho loci; skip geneid rows for non-ortho genes.
    for annotation_key, local_id in conn.execute(
        """
        SELECT annotation_key, local_id
        FROM gene_hit
        WHERE namespace = 'geneid'
        """
    ):
        geneid_rows_scanned += 1
        key = (int(annotation_key), int(local_id))
        if key in ortho_loci:
            loci_with_geneid.add(key)
        if geneid_rows_scanned % 5_000_000 == 0:
            print(
                f"# loci pass2: {geneid_rows_scanned:,} geneid rows, "
                f"{len(loci_with_geneid):,} ortho loci with geneid…",
                file=sys.stderr,
                flush=True,
            )

    print(
        f"# loci pass2 done: {geneid_rows_scanned:,} geneid rows scanned "
        f"({time.monotonic() - t1:.1f}s)",
        file=sys.stderr,
        flush=True,
    )

    accessions_with_any_geneid_locus: set[str] = set()
    accessions_never: set[str] = set(distinct_accessions)
    loci_total = len(ortho_loci)
    loci_hit = 0
    for locus, accs in ortho_loci.items():
        if locus in loci_with_geneid:
            loci_hit += 1
            accessions_with_any_geneid_locus.update(accs)
            accessions_never.difference_update(accs)

    never_sorted = sorted(accessions_never)
    return {
        "placement_rows": placement_rows,
        "loci_total": loci_total,
        "loci_with_geneid": loci_hit,
        "loci_missing_geneid": loci_total - loci_hit,
        "accessions_in_gene_hit": len(distinct_accessions),
        "accessions_with_geneid_locus": len(accessions_with_any_geneid_locus),
        "accessions_never_with_geneid_locus": len(never_sorted),
        "sample_never": never_sorted[:sample],
        "geneid_rows_scanned": geneid_rows_scanned,
    }


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
        "--sample",
        type=int,
        default=15,
        help="max missing accessions to print (default: 15)",
    )
    parser.add_argument(
        "--skip-loci",
        action="store_true",
        help="skip gene_hit locus co-occurrence (xref_meta only)",
    )
    args = parser.parse_args(argv)

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    try:
        # --- meta: accession equality ---
        t0 = time.monotonic()
        total, with_eq, missing_sample = audit_meta(conn, sample=args.sample)
        missing_eq = total - with_eq
        print("=== xref_meta: ncbiortholog accession == geneid accession ===")
        print(f"ncbiortholog accessions:     {total:,}")
        print(
            f"also present as geneid:      {with_eq:,}  ({_pct(with_eq, total)})"
        )
        print(
            f"missing equal geneid:        {missing_eq:,}  "
            f"({_pct(missing_eq, total)})"
        )
        print(f"query time:                  {time.monotonic() - t0:.1f}s")
        _print_sample(
            "sample missing equal geneid",
            missing_sample,
            args.sample,
        )
        print()

        if args.skip_loci:
            print("# skipped gene_hit locus audit (--skip-loci)")
            return 0

        # --- gene_hit: same locus co-occurrence ---
        print("=== gene_hit: same-locus co-occurrence with any geneid ===")
        stats = audit_loci(conn, sample=args.sample)
        loci_total = int(stats["loci_total"])
        loci_hit = int(stats["loci_with_geneid"])
        loci_miss = int(stats["loci_missing_geneid"])
        acc_total = int(stats["accessions_in_gene_hit"])
        acc_hit = int(stats["accessions_with_geneid_locus"])
        acc_miss = int(stats["accessions_never_with_geneid_locus"])

        print(f"ncbiortholog gene_hit rows:  {stats['placement_rows']:,}")
        print(f"distinct ortho loci:         {loci_total:,}")
        print(
            f"loci also having geneid:     {loci_hit:,}  "
            f"({_pct(loci_hit, loci_total)})"
        )
        print(
            f"loci missing geneid:         {loci_miss:,}  "
            f"({_pct(loci_miss, loci_total)})"
        )
        print(f"distinct ortho accessions:   {acc_total:,}")
        print(
            f"accessions with ≥1 geneid locus: {acc_hit:,}  "
            f"({_pct(acc_hit, acc_total)})"
        )
        print(
            f"accessions never with geneid:    {acc_miss:,}  "
            f"({_pct(acc_miss, acc_total)})"
        )
        print(f"geneid gene_hit rows scanned:  {stats['geneid_rows_scanned']:,}")
        _print_sample(
            "sample accessions never co-occurring with geneid",
            stats["sample_never"],  # type: ignore[arg-type]
            args.sample,
        )
    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
