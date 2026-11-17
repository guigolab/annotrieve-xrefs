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
ap.add_argument("--gate-check", type=int, default=40, help="shards per group (0 = skip)")
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


# ---------------------------------------------------------------------------
# Gate check: did the symbol gate reject names that are not locus tags?
# ---------------------------------------------------------------------------
import sys
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gene_corpus.stream.symbol_gate import is_locus_shaped, symbol_ok  # noqa: E402


def gate_check(shard: Path):
    """Return (n_genes_no_symbol, n_suspicious, examples) for one shard, streaming."""
    c = sqlite3.connect(f"file:{shard}?mode=ro", uri=True)
    cur = c.execute(
        "SELECT g.source_gene_id, g.primary_name, "
        "(SELECT x.accession FROM gene_xref x WHERE x.local_id=g.local_id AND x.namespace='locus_tag' LIMIT 1) "
        "FROM gene g WHERE g.primary_name IS NOT NULL AND NOT EXISTS "
        "(SELECT 1 FROM gene_xref s WHERE s.local_id=g.local_id AND s.namespace='symbol')")
    n = sus = 0
    ex: list[str] = []
    for sid, pn, lt in cur:
        n += 1
        p = pn.strip()
        echo = (lt and p == str(lt).strip()) or (sid and p == str(sid).strip()) \
            or (sid and p == str(sid).split("-", 1)[-1])
        # rejected only because of the locus-shaped patterns, not an ID/locus_tag echo
        if not echo and not symbol_ok(p) and is_locus_shaped(p):
            sus += 1
            if len(ex) < 6:
                ex.append(p)
    c.close()
    return n, sus, ex


b = a
if b.gate_check:
    print("\n=== gate check: genes without a symbol whose primary_name was rejected by the "
          "locus-shape patterns and is NOT an echo of locus_tag/ID ===")
    for label, cond in (("NAMED shards", "n_name_rows>0"), ("NAMELESS shards", "n_name_rows=0")):
        rows = m.execute("SELECT annotation_id, source_database FROM annotation WHERE "
                         f"status='indexed' AND in_catalogue=1 AND {cond}").fetchall()
        random.Random(a.seed).shuffle(rows)
        tot_n = tot_s = 0
        shards_with = 0
        by_src: Counter = Counter()
        seen: list[tuple[str, str, list[str]]] = []
        for aid, src in rows[:b.gate_check]:
            try:
                n, s, ex = gate_check(a.work_dir / "per_gff" / aid / "genes.sqlite")
            except Exception as e:
                print(f"  {aid}: ERROR {e}"); continue
            tot_n += n; tot_s += s
            if s:
                shards_with += 1
                by_src[src] += s
                if len(seen) < 12:
                    seen.append((aid, src, ex))
        print(f"\n{label}: {min(b.gate_check, len(rows))} shards, genes without symbol={tot_n:,}, "
              f"suspicious={tot_s:,}, shards_with_suspicious={shards_with}, by_source={dict(by_src)}")
        for aid, src, ex in seen:
            print(f"  {aid} {src} {ex}")
