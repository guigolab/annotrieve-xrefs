"""
Smoke test of /sequences on the real index and GFF files (read-only, no API server).

    python3 tools/smoke_sequences.py \
        --serve-dir /data/annotrieve/gene_corpus/serve --annotations-dir /data/annotrieve/files

Checks, with timings and peak RSS (all local; no NCBI and no gffread unless --export):
  1. count of a query (genes, base-pair upper bound);
  2. locus_ref lookups through loci_ref.sqlite: latency of N random refs, all must resolve;
  3. resolving N genes against the real GFF files: skip reasons, candidates per gene, and the
     transcript chosen by several presets (shows how often the presets disagree);
  4. optional --export N: a real export of N genes through gffread and NCBI into a temp folder.
"""
import argparse
import random
import resource
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "api"), str(ROOT)]
import settings  # noqa: E402

PRESETS = ("longest", "longest_spliced", "longest_cds", "shortest_cds", "most_exons", "fewest_exons")


def rss() -> str:
    """Current (not peak) resident memory of this process, MiB."""
    with open("/proc/self/statm") as fh:
        pages = int(fh.read().split()[1])
    return f"rss={pages * 4096 / 1048576:.0f}MiB"


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    print(f"[{(time.perf_counter() - t) * 1000:9.1f} ms] {label}  ({rss()})")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--serve-dir", type=Path, required=True)
    ap.add_argument("--annotations-dir", type=Path, required=True)
    ap.add_argument("--names", default="tp53", help="comma-separated exact gene names")
    ap.add_argument("--taxid", type=int, default=0, help="restrict to this taxon (0 = all)")
    ap.add_argument("--genes", type=int, default=40, help="genes to resolve against the GFFs")
    ap.add_argument("--refs", type=int, default=1000, help="random locus_refs to look up")
    ap.add_argument("--export", type=int, default=0,
                    help="also export this many genes through gffread and NCBI (needs network)")
    ap.add_argument("--profile", action="store_true",
                    help="profile the GFF resolving loop (top functions by cumulative time)")
    ap.add_argument("--action", default="cds", choices=["cds", "transcript", "protein"])
    a = ap.parse_args(argv)

    s = settings.settings
    s.PRESENCE_DIR, s.ANNOTATIONS_DIR = str(a.serve_dir), str(a.annotations_dir)

    from helpers.transcript_select import parse_rule, select_transcripts
    from services.seq_request import LociEntry, SeqRequest, parse_request
    from services.seq_resolve import GeneResolver
    from services.seq_source import Source

    query = {"names": a.names.split(","), "alias": "off"}
    if a.taxid:
        query["include_taxids"] = [a.taxid]
    req = parse_request({"query": query, "action": a.action, "max_genes": 100_000})
    src = timed("open index connections + build query", lambda: Source(req))
    print("data_version", src.data_version, "| terms", src.terms())

    counted = timed("count genes", lambda: src.count(100_000))
    print(f"  genes={counted.n_genes:,} complete={counted.complete} "
          f"bp_upper_bound={counted.bp_upper_bound:,}")
    genes = [g for _r, g, _x in src.items() if g is not None]
    if not genes:
        print("no genes matched; nothing more to check")
        return 1
    random.Random(1).shuffle(genes)

    print(f"\n--- locus_ref lookup ({min(a.refs, len(genes))} refs) ---")
    entries = tuple(LociEntry(g.annotation_id, g.local_id) for g in genes[: a.refs])
    ref_req = SeqRequest(req.action, req.rule, "longest", 0, False, "skip", 100_000, None, entries)
    ref_src = Source(ref_req)
    looked = timed(f"look up {len(entries)} locus_refs", lambda: list(ref_src.items()))
    miss = [r for r, g, why in looked if g is None]
    by_ref = {g.locus_ref: g for g in genes}
    wrong = [r for r, g, _ in looked if g is not None and
             (g.seqid, g.start, g.end, g.gene_id) != (by_ref[r].seqid, by_ref[r].start,
                                                      by_ref[r].end, by_ref[r].gene_id)]
    print(f"  unresolved={len(miss)} different_from_query_row={len(wrong)}")
    ref_src.close()

    n = min(a.genes, len(genes))
    print(f"\n--- resolving {n} genes against the real GFF files (action={a.action}) ---")
    resolver = GeneResolver(a.action, req.rule, 0)
    reasons: Counter = Counter()
    n_cands: Counter = Counter()
    disagree: Counter = Counter()
    shown = 0
    prof = None
    if a.profile:
        import cProfile
        prof = cProfile.Profile()
        prof.enable()
    t0 = time.perf_counter()
    for g in genes[:n]:
        res = resolver.resolve(g)
        reasons[res.reason or "ok"] += 1
        n_cands[min(res.n_candidates, 10)] += 1
        if res.candidates:
            picks = {p: tuple(c.transcript_id for c in select_transcripts(
                res.candidates, parse_rule(p, a.action), a.action).chosen) for p in PRESETS}
            disagree[len(set(picks.values()))] += 1
            if shown < 3:
                shown += 1
                print(f"  {g.locus_ref} {g.symbol} {g.seqid}:{g.start}-{g.end} "
                      f"candidates={res.n_candidates} reason={res.reason or 'ok'}")
                for c in res.candidates[:4]:
                    print(f"     {c.transcript_id:<28} {c.metrics()}")
                print(f"     presets -> {picks}")
        elif shown < 6 and res.reason:
            shown += 1
            print(f"  {g.locus_ref} {g.symbol}: skipped {res.reason}: {res.message}")
    if prof is not None:
        prof.disable()
    print(f"[{(time.perf_counter() - t0) * 1000:9.1f} ms] {n} genes resolved "
          f"({(time.perf_counter() - t0) * 1000 / n:.1f} ms/gene)  ({rss()})")
    if prof is not None:
        import io
        import pstats
        buf = io.StringIO()
        pstats.Stats(prof, stream=buf).sort_stats("cumulative").print_stats(14)
        print("\n".join(ln for ln in buf.getvalue().splitlines() if ln.strip())[:3500])
    print(f"  open tabix handles before close: {len(resolver.tabix)}")
    resolver.close()
    print("  outcome:", dict(reasons))
    print("  candidate transcripts per gene (10 = 10 or more):", dict(sorted(n_cands.items())))
    print("  number of different transcripts chosen by the 6 presets:", dict(sorted(disagree.items())))
    src.close()

    if a.export:
        from services.seq_export import run_export
        print(f"\n--- real export of {a.export} genes through gffread and NCBI ---")
        exp_req = parse_request({"query": query, "action": a.action, "max_genes": 100_000,
                                 "on_error": "skip"})
        exp_src = Source(exp_req)
        out = Path(tempfile.mkdtemp(prefix="smoke-seq-"))

        def go():
            items = (it for i, it in enumerate(exp_src.items()) if i < a.export)
            return run_export(exp_req, items, out_dir=out, total=a.export)

        result = timed("export", go)
        print(f"  sequences={result.n_sequences} genes={result.n_genes} skipped={result.n_skipped} "
              f"reasons={result.reasons}\n  files in {out}")
        for line in result.fasta_path.read_text().splitlines()[:4]:
            print("  ", line[:140])
        exp_src.close()

    print(f"\npeak RSS MiB {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
