"""
Why do some indexed shards have no searchable names?  (read-only, bounded memory)

    python3 tools/diagnose_nameless.py --work-dir /data/annotrieve/gene_corpus --sample 40

1. Breakdown of named vs nameless indexed shards by source_database (manifest only).
2. For a random sample of nameless shards: genes, genes with primary_name, xref rows per
   namespace, and up to 8 example primary_names, so a human can judge whether the gate
   (stream/symbol_gate.py) dropped real symbols.
"""
import argparse, random, sqlite3
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--work-dir", type=Path, required=True)
ap.add_argument("--sample", type=int, default=40)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

m = sqlite3.connect(f"file:{a.work_dir/'serve'/'manifest.sqlite'}?mode=ro", uri=True)
print("source_database  indexed  nameless  pct_nameless  (in_catalogue only)")
for src, n, z in m.execute(
    "SELECT source_database, COUNT(*), SUM(n_name_rows=0) FROM annotation "
    "WHERE status='indexed' AND in_catalogue=1 GROUP BY source_database"):
    print(f"{src:<16} {n:>7} {z:>9} {100*z/n:>9.1f}%")
print("\nnameless by n_genes bucket:")
for lab, n in m.execute(
    "SELECT CASE WHEN n_genes<100 THEN '<100' WHEN n_genes<10000 THEN '100-10k' "
    "ELSE '>=10k' END, COUNT(*) FROM annotation WHERE status='indexed' AND in_catalogue=1 "
    "AND n_name_rows=0 GROUP BY 1"):
    print(f"  {lab:<8} {n}")

ids = [r for r in m.execute(
    "SELECT annotation_id, source_database, organism_name, n_genes FROM annotation "
    "WHERE status='indexed' AND in_catalogue=1 AND n_name_rows=0")]
random.Random(a.seed).shuffle(ids)
print(f"\nsample of {min(a.sample, len(ids))} nameless shards:")
for aid, src, org, ng in ids[:a.sample]:
    p = a.work_dir / "per_gff" / aid / "genes.sqlite"
    try:
        c = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
        genes = c.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
        with_pn = c.execute("SELECT COUNT(*) FROM gene WHERE primary_name IS NOT NULL").fetchone()[0]
        ns = dict(c.execute("SELECT namespace, COUNT(*) FROM gene_xref GROUP BY namespace ORDER BY 2 DESC LIMIT 6"))
        ex = [r[0] for r in c.execute("SELECT primary_name FROM gene WHERE primary_name IS NOT NULL LIMIT 8")]
        c.close()
    except Exception as e:
        print(f"{aid} {src}: ERROR {e}"); continue
    print(f"{aid} {src} {org!r} genes={genes} with_primary_name={with_pn} xref_ns={ns} examples={ex}")
