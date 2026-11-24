"""
Worst-case memory check of the search API (API-LOCKED-DECISIONS.md, sections 3 and 4).

    python3 tools/check_api_memory.py --scratch /data/annotrieve-xrefs/sequences/tmp [--limit-mb 486]

Builds a synthetic presence index with 18,000 annotations in which 50 gene names are present in
every annotation (the heaviest possible /hits answer), starts the real app in this process, then
replays the worst requests from many client threads at once: /hits with 50 names in both shapes
and as TSV, /names, /health. The API gates decide how many run together, exactly as in production.
Exits 1 when the peak RSS of this process (one worker) is above --limit-mb.

Memory of this run: the synthetic index is about 3 MB; the process peak is what is reported.
The scratch directory must be on the network drive (never local disk); it is removed at the end.
"""
import argparse
import concurrent.futures as cf
import resource
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "api"), str(ROOT)]

N_ANN = 18_000
N_NAMES = 50


def rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # Linux: KiB


def build_index(path: Path) -> list[str]:
    from gene_corpus.sql.postings import encode_postings
    from gene_corpus.sql.serve_presence import _OUT_SCHEMA

    conn = sqlite3.connect(path)
    conn.executescript(_OUT_SCHEMA)
    conn.executemany(
        "INSERT INTO annotation (ord, annotation_id, taxid, organism_name, assembly_accession, "
        "assembly_name, source, release_date, gff_path) VALUES (?,?,?,?,?,?,?,?,?)",
        [(o, f"{o:032x}", o % 9000 + 1, f"Genus species {o}", f"GCF_{o:09d}.1", f"asm{o}",
          1 + o % 3, "2024-01-01", "/x/y/z.gff.gz") for o in range(1, N_ANN + 1)])
    conn.executemany("INSERT OR IGNORE INTO lineage VALUES (?,?)",
                     [(t, o) for o in range(1, N_ANN + 1) for t in (o % 9000 + 1, 40674, 2759)])
    everyone = list(range(1, N_ANN + 1))
    blob = encode_postings(everyone, N_ANN)
    names = [f"gene{i:03d}" for i in range(N_NAMES)]
    for i, name in enumerate(names, 1):
        conn.execute("INSERT INTO name_dict VALUES (?,?,?,?)", (name, i, N_ANN, 0))
        conn.execute("INSERT INTO posting VALUES (?,?,?)", (i, blob, None))
    conn.executemany("INSERT INTO meta VALUES (?,?)",
                     [("data_version", "memcheck"), ("schema_version", "1")])
    conn.commit()
    conn.close()
    return names


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", type=Path, required=True, help="directory on the network drive")
    ap.add_argument("--limit-mb", type=float, default=486.0,
                    help="peak RSS allowed for one worker (budget table, section 4)")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--clients", type=int, default=24)
    a = ap.parse_args()

    a.scratch.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="memcheck-", dir=a.scratch))
    try:
        serve = work / "serve"
        (serve / "v1").mkdir(parents=True)
        names = build_index(serve / "v1" / "presence.sqlite")
        (serve / "current").symlink_to("v1")

        import settings as settings_mod
        settings_mod.settings.PRESENCE_DIR = str(serve)
        from fastapi.testclient import TestClient
        from main import create_app

        base = rss_mb()
        urls = []
        for shape in ("compact", "records"):
            urls.append(("/hits", {"names": ",".join(names), "shape": shape}))
        urls.append(("/hits", {"names": ",".join(names), "format": "tsv"}))
        urls.append(("/hits", {"names": names[0], "shape": "records"}))
        urls.append(("/names", {"prefix": "gene"}))
        urls.append(("/health", {}))

        codes: dict[int, int] = {}
        t0 = time.monotonic()
        with TestClient(create_app()) as client:
            def one(i):
                path, params = urls[i % len(urls)]
                r = client.get(path, params=params)
                return r.status_code, len(r.content)

            with cf.ThreadPoolExecutor(max_workers=a.clients) as pool:
                for status, _size in pool.map(one, range(a.rounds * a.clients)):
                    codes[status] = codes.get(status, 0) + 1
        peak = rss_mb()
        print(f"requests: {sum(codes.values())} in {time.monotonic() - t0:.1f} s, "
              f"status counts {dict(sorted(codes.items()))}")
        print(f"RSS before replay (app loaded): {base:.0f} MB; peak: {peak:.0f} MB; "
              f"limit for one worker: {a.limit_mb:.0f} MB")
        bad = {c: n for c, n in codes.items() if c not in (200, 503)}
        if bad:
            print(f"FAIL: unexpected status codes {bad}")
            return 1
        if codes.get(200, 0) == 0:
            print("FAIL: no request succeeded")
            return 1
        if peak > a.limit_mb:
            print("FAIL: peak above the limit")
            return 1
        print("OK")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
