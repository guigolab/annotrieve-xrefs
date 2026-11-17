"""
Dry run: which shards would gain symbols under the relaxed symbol gate?  (read-only)

    python3 tools/measure_gate_relax.py --work-dir /data/annotrieve/gene_corpus \
        [--sample 300] [--workers 2] [--out gate_gain.tsv]

For each indexed, in-catalogue shard, looks at genes with a primary_name but no `symbol`
xref row, drops names that merely echo the gene's ID / locus_tag (including the Ensembl
`gene:` prefix), and counts those the new `symbol_ok` accepts. Only symbols can be
measured here (primary_name is Name-first); rejected aliases are not stored in shards, so
alias effects need a GFF re-scan and are not measured.

Output TSV (`--out`): annotation_id, source_database, gained_genes, examples. The ids
with gained_genes > 0 are the re-harvest candidates (feed them to `--ids`).
Memory: streamed cursor, a handful of counters per shard; peak RSS is reported.
"""
import argparse, random, resource, sqlite3, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gene_corpus.stream.symbol_gate import symbol_ok  # noqa: E402


def _echo(p: str, sid: str | None, lt: str | None) -> bool:
    if lt and p == str(lt).strip():
        return True
    if sid:
        s = str(sid).strip()
        cands = {s, s.split(":", 1)[-1], s[5:] if s.startswith("gene-") else s}
        if p in cands:
            return True
    return False


def scan(task):
    work, aid, src = task
    shard = Path(work) / "per_gff" / aid / "genes.sqlite"
    gained = 0
    ex = []
    try:
        c = sqlite3.connect(f"file:{shard}?mode=ro", uri=True)
        cur = c.execute(
            "SELECT g.source_gene_id, g.primary_name, "
            "(SELECT x.accession FROM gene_xref x WHERE x.local_id=g.local_id AND x.namespace='locus_tag' LIMIT 1) "
            "FROM gene g WHERE g.primary_name IS NOT NULL AND NOT EXISTS "
            "(SELECT 1 FROM gene_xref s WHERE s.local_id=g.local_id AND s.namespace='symbol')")
        for sid, pn, lt in cur:
            p = pn.strip()
            if not p or _echo(p, sid, lt):
                continue
            if symbol_ok(p):
                gained += 1
                if len(ex) < 5:
                    ex.append(p)
        c.close()
    except Exception as e:  # report, never stop the run
        return aid, src, -1, [f"ERROR {type(e).__name__}: {e}"]
    return aid, src, gained, ex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--sample", type=int, default=0, help="random N shards (0 = all)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    m = sqlite3.connect(f"file:{a.work_dir/'serve'/'manifest.sqlite'}?mode=ro", uri=True)
    rows = m.execute("SELECT annotation_id, source_database FROM annotation "
                     "WHERE status='indexed' AND in_catalogue=1 ORDER BY ord").fetchall()
    if a.sample:
        rows = random.Random(a.seed).sample(rows, min(a.sample, len(rows)))
    tasks = [(str(a.work_dir), aid, src) for aid, src in rows]
    t0 = time.time()
    res = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(scan, tasks, chunksize=2), 1):
            res.append(r)
            if i % 200 == 0:
                print(f"  {i}/{len(tasks)} {i/(time.time()-t0):.1f} shards/s", file=sys.stderr, flush=True)
    errs = [r for r in res if r[2] < 0]
    hit = [r for r in res if r[2] > 0]
    tot = sum(r[2] for r in hit)
    print(f"shards scanned={len(res):,} errors={len(errs)} shards_with_gain={len(hit):,} "
          f"({100*len(hit)/max(1,len(res)):.1f}%) genes_gained={tot:,} "
          f"time={time.time()-t0:.0f}s parent_peak_rss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f}MiB")
    by = {}
    for _aid, src, g, _ in hit:
        d = by.setdefault(src, [0, 0]); d[0] += 1; d[1] += g
    for src, (n, g) in sorted(by.items()):
        print(f"  {src:<10} shards={n:<6} genes={g:,}")
    for aid, src, g, ex in sorted(hit, key=lambda r: -r[2])[:15]:
        print(f"  {aid} {src} +{g:,} {ex}")
    for r in errs[:5]:
        print("  ERR", r)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("annotation_id\tsource_database\tgained_genes\texamples\n")
            for aid, src, g, ex in hit:
                fh.write(f"{aid}\t{src}\t{g}\t{','.join(ex)}\n")
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
