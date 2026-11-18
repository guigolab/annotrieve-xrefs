"""
Smoke test of the real presence index through the API service code (read-only).

    python3 tools/smoke_presence.py --serve-dir /data/annotrieve/gene_corpus/serve

Opens <serve-dir>/current/presence.sqlite exactly as the API does and prints answers and
timings for a few queries, plus a sample of alias-only names. No server, no old tables.
"""
import argparse, random, resource, sqlite3, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "api"), str(ROOT)]
import settings  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--serve-dir", type=Path, required=True)
a = ap.parse_args()
settings.settings.PRESENCE_DIR = str(a.serve_dir)

from helpers.presence import presence_reader  # noqa: E402
from services.presence_hits import list_name_hits_presence, list_names_by_prefix_presence  # noqa: E402


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    print(f"[{(time.perf_counter() - t) * 1000:8.1f} ms] {label}")
    return out


r = timed("open reader (builds annotation mask)", presence_reader)
print("data_version", r.data_version, "annotations in mask", r.all_mask.bit_count())

for names in ("tp53", "tp53,brca1,brca2", "selenot,selenop,selenow", "polq", "e3p82_00001"):
    b = timed(f"/hits names={names}", lambda n=names: list_name_hits_presence(r, n))
    print("   total", b["total"], [(t["name"], t["n_hits"], t["n_alias_extra"]) for t in b["terms"]])

b = timed("/hits tp53 include_taxids=40674 (Mammalia)",
          lambda: list_name_hits_presence(r, "tp53", include_taxids="40674"))
print("   total", b["total"])
b = timed("/hits tp53 alias=include", lambda: list_name_hits_presence(r, "tp53", alias="include"))
print("   total", b["total"])

for prefix in ("seleno", "tp53", "brca"):
    try:
        n = timed(f"/names prefix={prefix}", lambda p=prefix: list_names_by_prefix_presence(r, p))
        print("   total", n["total"], [(x["name"], x["n_annotations"], x["n_alias_extra"]) for x in n["results"][:5]])
    except Exception as e:
        print("   ", type(e).__name__, getattr(e, "code", ""), getattr(e, "extra", ""))

n = timed("/names prefix=seleno scoped to 9606", lambda: list_names_by_prefix_presence(r, "seleno", include_taxids="9606"))
print("   total", n["total"])

c = sqlite3.connect(f"file:{r.path}?mode=ro", uri=True)
rows = c.execute("SELECT name FROM name_dict WHERE n_sym=0 AND name_id % 5000 = 0 LIMIT 25").fetchall()
print("sample of alias-only names:", [x[0] for x in rows])
print("peak RSS MiB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024)


# ---------------------------------------------------------------------------
# /loci on the real loci.sqlite (cold NFS read, then warm)
# ---------------------------------------------------------------------------
from helpers.presence import annotation_info, open_request_connections  # noqa: E402
from services.presence_loci import build_query, estimate, iter_loci  # noqa: E402
from services.loci_formats import stream_tsv  # noqa: E402

print("\n--- /loci ---")
pres, loci, dv = open_request_connections()
INFO = annotation_info(r)
print("loci data_version", dv, "rows", loci.execute("SELECT value FROM meta WHERE key='n_loci'").fetchone()[0])


def loci_query(names, **kw):
    q = build_query(pres, r.all_mask, dv, names, **kw)
    est = estimate(pres, loci, q)
    t = time.perf_counter()
    first = None
    n = 0
    for row in iter_loci(pres, loci, q, None, INFO):
        if first is None:
            first = (time.perf_counter() - t) * 1000
        n += 1
    total = (time.perf_counter() - t) * 1000
    print(f"[{total:9.1f} ms, first row {first or 0:7.1f} ms] /loci names={names} {kw or ''} -> "
          f"{n:,} rows (annotations {est['n_annotations']:,}, upper bound {est['n_loci_upper_bound']:,})")
    return q


for names in ("tp53", "tp53", "brca2", "polq", "tp53,brca1,brca2,polq,actb"):
    loci_query(names)
loci_query("tp53", include_taxids="40674")
loci_query("tp53", alias="off")
q = loci_query("seleno1,selenop,selenow")

t = time.perf_counter()
size = sum(len(chunk) for chunk in stream_tsv(q, iter_loci(pres, loci, q, None, INFO), 2_000_000))
print(f"[{(time.perf_counter() - t) * 1000:9.1f} ms] tsv stream of the last query: {size / 1024:.0f} KiB")
row = next(iter_loci(pres, loci, loci_query("tp53", alias="off"), None, INFO))
print("example row:", {k: row[k] for k in ("annotation_id", "symbol", "seqid", "start", "end", "strand", "gene_id", "gff_url")})
print("peak RSS MiB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024)
