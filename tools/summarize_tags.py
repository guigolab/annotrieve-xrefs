"""
Summarise the per-file TSV written by measure_gff_tags.py --out, with organism names.

    python3 tools/summarize_tags.py /tmp/tags_all.tsv \
        --manifest /data/annotrieve/gene_corpus/serve/manifest.sqlite

The scanned file name (`RefSeq_<hash>.gff.gz`) is joined to the manifest `catalogue` table by the
basename of `bgzip_path`. The script prints how many files matched, so a wrong join is visible.
Read-only; memory is one row per scanned file.
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", type=Path)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--show", type=int, default=40, help="max species listed per group")
    a = ap.parse_args()

    conn = sqlite3.connect(f"file:{a.manifest}?mode=ro", uri=True)
    info: dict[str, tuple[str, int]] = {}
    for path, organism, taxid in conn.execute(
            "SELECT bgzip_path, organism_name, taxid FROM catalogue WHERE bgzip_path <> ''"):
        info[path.rsplit("/", 1)[-1]] = (organism or "?", int(taxid or 0))

    rows = []
    with open(a.tsv) as fh:
        next(fh)
        for line in fh:
            name, cls, genes, tx, flagged = line.rstrip("\n").split("\t")
            organism, taxid = info.get(name, ("(not in catalogue)", 0))
            rows.append((name, cls, organism, taxid, int(genes), int(tx), int(flagged)))
    if not rows:
        sys.exit("empty TSV")

    matched = sum(1 for r in rows if r[2] != "(not in catalogue)")
    print(f"{len(rows)} files, {matched} matched to the catalogue by file name")

    by_class: dict[str, dict] = defaultdict(lambda: {"files": 0, "flagged_files": 0,
                                                     "species": set(), "flagged_species": set()})
    for _n, cls, organism, _t, _g, _tx, flagged in rows:
        d = by_class[cls]
        d["files"] += 1
        d["species"].add(organism)
        if flagged:
            d["flagged_files"] += 1
            d["flagged_species"].add(organism)
    print(f"\n{'source':<12}{'files':>8}{'flagged':>9}{'species':>9}{'flagged species':>17}")
    for cls, d in sorted(by_class.items()):
        print(f"{cls:<12}{d['files']:>8}{d['flagged_files']:>9}{len(d['species']):>9}"
              f"{len(d['flagged_species']):>17}")

    for cls in sorted(by_class):
        per_species: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
        for _n, c, organism, _t, genes, _tx, flagged in rows:
            if c == cls and flagged:
                s = per_species[organism]
                s[0] += 1
                s[1] += genes
                s[2] += flagged
        if per_species:
            print(f"\n-- flagged species in {cls} (files, genes, flagged genes) --")
            for organism, (f, g, fl) in sorted(per_species.items())[: a.show]:
                print(f"  {organism:<45}{f:>4} files {g:>9,} genes {fl:>9,} flagged ({100 * fl // max(g, 1)}%)")


if __name__ == "__main__":
    main()
