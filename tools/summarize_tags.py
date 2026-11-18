"""
Summarise the per-file TSV written by measure_gff_tags.py --out, with organism names.

    python3 tools/summarize_tags.py /tmp/tags_refseq.tsv \
        --presence /data/annotrieve/gene_corpus/serve/current/presence.sqlite

File names look like `RefSeq_<hash>.gff.gz`; the hash is assumed to be the annotation_id
(the script prints how many names matched, so a wrong assumption is visible).
Read-only; memory is one row per scanned file.
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", type=Path)
    ap.add_argument("--presence", type=Path, required=True)
    ap.add_argument("--show", type=int, default=60, help="max organisms listed per group")
    a = ap.parse_args()

    conn = sqlite3.connect(f"file:{a.presence}?mode=ro", uri=True)
    org = dict(conn.execute("SELECT annotation_id, organism_name FROM annotation"))

    rows = []
    with open(a.tsv) as fh:
        next(fh)
        for line in fh:
            name, cls, genes, tx, flagged = line.rstrip("\n").split("\t")
            h = name.split("_", 1)[1].removesuffix(".gz").removesuffix(".gff")
            rows.append((h, cls, int(genes), int(tx), int(flagged)))
    if not rows:
        sys.exit("empty TSV")

    matched = sum(1 for r in rows if r[0] in org)
    print(f"{len(rows)} files, {matched} matched to an organism by annotation_id")
    flagged = [r for r in rows if r[4] > 0]
    print(f"files with a canonical flag: {len(flagged)} / {len(rows)} "
          f"({100 * len(flagged) / len(rows):.1f}%)")
    tot_genes = sum(r[2] for r in rows)
    print(f"genes in scanned files {tot_genes:,}; genes with a flagged transcript "
          f"{sum(r[4] for r in rows):,}")

    def fmt(r):
        pct = 100 * r[4] / r[2] if r[2] else 0
        return f"  {org.get(r[0], '(not in index)'):<45} genes={r[2]:>7,} flagged={r[4]:>7,} ({pct:.0f}%)"

    print("\n-- FLAGGED --")
    for r in sorted(flagged, key=lambda r: org.get(r[0], ""))[: a.show]:
        print(fmt(r))
    print("\n-- NOT FLAGGED (first few) --")
    for r in sorted([r for r in rows if r[4] == 0], key=lambda r: org.get(r[0], ""))[: a.show]:
        print(fmt(r))


if __name__ == "__main__":
    main()
