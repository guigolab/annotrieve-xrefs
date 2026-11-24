"""Load replay against a running API container.

    python3 tools/load_check_container.py <seconds> <clients>   # API on 127.0.0.1:15000, index from check_api_memory.build_index
"""
import sys, time, threading, collections, httpx
base = "http://127.0.0.1:15000"
names = ",".join(f"gene{i:03d}" for i in range(50))
reqs = [("/hits", {"names": names, "shape": "records"}), ("/hits", {"names": names}), ("/hits", {"names": names, "format": "tsv"}),
        ("/hits", {"names": "gene000", "shape": "records"}), ("/names", {"prefix": "gene"}), ("/names", {"prefix": "zz0000"}), ("/health", {})]
dur = float(sys.argv[1]); clients = int(sys.argv[2])
codes = collections.Counter(); lat = collections.defaultdict(list); stop = time.time() + dur; lock = threading.Lock()
def work(k):
    c = httpx.Client(timeout=60); i = k
    while time.time() < stop:
        path, params = reqs[i % len(reqs)]; i += 1
        params = dict(params)
        if path == "/hits" and i % 3 == 0: params["include_taxids"] = str(40674 if i % 2 else 2759)  # vary URL: no nginx here, but vary scope
        t = time.perf_counter(); r = c.get(base + path, params=params); d = (time.perf_counter() - t) * 1000
        with lock: codes[r.status_code] += 1; lat[path + ("?" + params.get("shape", params.get("format", ""))) ].append(d)
ts = [threading.Thread(target=work, args=(k,)) for k in range(clients)]
[t.start() for t in ts]; [t.join() for t in ts]
print("status", dict(codes))
for k, v in sorted(lat.items()):
    v.sort(); print(f"{k:28s} n={len(v):5d} p50={v[len(v)//2]:7.1f} p95={v[int(len(v)*.95)]:7.1f} max={v[-1]:7.1f} ms")
