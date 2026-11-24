# API: locked decisions

Status: decisions locked 2026-11-24. This file tracks the work: tick the boxes in section 7 and add a line to the progress log (section 8) with each change.

Scope: only the search API (`/names`, `/hits`, `/loci`, `/annotations/{id}/...`) and the files it reads. Sequence export logic (`/sequences`) is out of scope except for its memory use, because it runs in the same container.

Goal: a research-grade API that stays fast under heavy use, never exceeds 2 GB of RAM, and writes almost nothing to local disk.

## 1. Hard rules

1. **RAM: at most 2 GB for the whole API container** (both workers, shared memory and child processes). This is enforced by the container limit, not just by good intentions (section 4).
2. **Local disk: none by design.** Everything that is more than a few KB goes to the NFS data volume already mounted in the compose file. Section 5 lists what still touches local disk and how each is removed.
3. **Python stays.** Measured: the Python worker baseline is 66 MiB, request compute is 1 to 30 ms, `/loci` streaming costs 9 microseconds per row. The slow part was waiting on the network drive, which Python does not cause.

## 2. Facts these decisions rest on

All measured or read from the code on 2026-11-24.

| Fact | Value | Source |
|---|---|---|
| `presence.sqlite` (the search index) | 310 MB | build report |
| `loci.sqlite` | 5,319 MB, 82.2M rows | `measure-out/report.md` |
| Annotations | about 18,000 | build report |
| Python worker after importing the app | 66 MiB | measured |
| Gunicorn workers | 2 (`api/Dockerfile`) | code |
| Threads per worker | 40 (framework default, not overridden) | `main.py`, library default |
| SQLite page cache per presence connection | 256 MiB (default, not overridden in compose) | `settings.py` |
| `/hits` worst case, one request, whole corpus matched | 21 MB (compact) to 63 MB (records, 50 names) | measured, section 3 |
| `/hits` concurrent limit per worker today | 16 | `helpers/concurrency.py` |
| `/loci` concurrent limit per worker today | 4, 16 MiB cache each | `settings.py` |
| Server time of an uncached `/hits` today | 100 to 125 ms for large results; 10 to 20 ms for names already in memory | live API probes |
| `/health` | returns a fixed `ok`, checks nothing | `routers/health.py` |
| Container memory limit in compose | none | `docker-compose.yml` |

What this says about today: 16 concurrent worst-case `/hits` per worker is 16 x 63 MB = about 1 GB per worker, 2 GB for the container, before anything else. So **today the 2 GB rule is not guaranteed**, and the caps below are needed.

## 3. Worst-case memory of one request (measured)

`/hits` for a name present in all 18,000 annotations, building the response the way FastAPI does (service + encoder + JSON bytes), peak Python allocation:

| Request | Peak | Response size |
|---|---|---|
| 1 name, compact | 21 MB | 2.0 MB |
| 1 name, records | 44 MB | 6.4 MB |
| 50 names, compact | 40 MB | 4.5 MB |
| 50 names, records | 63 MB | 8.9 MB |

After the pre-built-response step (decision D4) the same requests need about the response size plus small overhead (a few MB).

## 4. Memory budget and how the 2 GB is guaranteed

### Three layers

1. **Hard limit in the compose file:** `mem_limit: 2g` and `memswap_limit: 2g` (no swap). The container can never use more. Memory in `/dev/shm` counts toward this limit.
2. **Caps on every consumer** (table below), so the sum of all worst cases stays under the limit.
3. **Automated check:** a script that replays the worst-case requests against a built index and fails if peak memory goes above the numbers in this table (`tools/check_api_memory.py`, to be written, see section 7).

If a bug ever pushes the container to its limit, the container is restarted by Docker; the host and the other services are not affected. Layers 2 and 3 are there so that this does not happen.

### Budget, worst case (everything at its cap at the same time)

| Item | Per worker | Both workers |
|---|---|---|
| Python baseline | 66 MB | 132 MB |
| Annotation list in memory (D3) | 7 MB | 14 MB |
| Presence connections: at most 12 threads, 2 MiB cache each (D3) | 24 MB | 48 MB |
| `/hits`: 3 concurrent x 63 MB worst case | 189 MB | 378 MB |
| `/loci`: 3 concurrent x (16 MiB cache + 10 MB buffers) | 78 MB | 156 MB |
| `/annotations/*`: 3 concurrent x (16 MiB cache + 10 MB) (D6) | 78 MB | 156 MB |
| `/sequences`, Python side: 32 MB cache + 1 job x 20 MB | 52 MB | 104 MB |
| **Workers subtotal** | **about 486 MB** | **about 972 MB** |
| Search index in `/dev/shm`, one copy (D1) | | 310 MB |
| Same, during a version switch (old + new copy, briefly) | | 620 MB |
| `gffread` child processes, 2 at most, 128 MB data limit each (D9) | | 256 MB |
| **Total at the worst moment** | | **about 1,850 MB** |

Margin at the worst moment: about 200 MB (10%). Typical use is far below this: the worst case needs the heaviest request type at its cap in every slot, plus a version switch, at once. After D4 (pre-built responses) the `/hits` line drops from 378 MB to about 60 MB and the total to about 1,530 MB.

### Settings that must change for the budget to hold

| Setting | Today | Locked value |
|---|---|---|
| `HITS_MAX_INFLIGHT` | 16 per worker | 3 per worker (raise after D4) |
| `LOCI_MAX_INFLIGHT` | 4 per worker | 3 per worker |
| `/annotations/*` concurrency gate | none | 3 per worker |
| `SEQUENCES_MAX_INFLIGHT` | 2 per worker | 1 per worker |
| `READ_CACHE_SIZE` (SQLite cache per presence connection) | 256 MiB | 2 MiB |
| SQLite `soft_heap_limit` | 256 MiB | 64 MiB |
| Connections | one per thread (40 threads per worker) | one per thread, thread pool cut to 12 per worker (`THREADPOOL_SIZE`) |
| Gunicorn `--max-requests` | not set | 5000 with jitter, so any slow leak is cleared |
| `SEQUENCES_MAX_INDEX_MB` | 0 (no limit) | 256 (opening a large index can use as much memory again) |
| Container `mem_limit` / `memswap_limit` | none | 2g / 2g |

## 5. Local disk: what is written, and where it goes

| Item | Today | Locked |
|---|---|---|
| SQLite temporary files (`temp_store=FILE`) | `/tmp`, a local Docker volume (`xrefs_tmp`) | `TMPDIR` and `SQLITE_TMPDIR` on NFS (`/data/results/tmp`); queries must also avoid sorts that spill (D7) |
| NCBI rate-limit lock file | `/tmp`, local | stays local: a few bytes, and file locking is unreliable on NFS |
| Container logs (`--access-logfile -`) | Docker default driver, unbounded, local | `logging: json-file` with `max-size: 10m`, `max-file: 3` (at most 30 MB) |
| Sequence results and scratch | NFS (`/data/results`) | unchanged |
| Search index copy | n/a | RAM (`/dev/shm`), not disk |
| Container image layers | local | unchanged, outside this scope |

## 6. Decisions

**D1. The whole search index lives in RAM.** At startup and whenever a new data version is published, copy `presence.sqlite` (310 MB) from NFS into `/dev/shm` once, shared by both workers (a lock stops them copying twice). All `/hits` and `/names` lookups then read memory, never the network drive. If the file is bigger than `PRESENCE_RAM_MAX_MB` (450) or `/dev/shm` is too small, the API keeps working from NFS and logs a warning. The compose file gets `shm_size: 1g` (room for two copies during a switch). Also: stop calling `stat` and symlink resolution on NFS for every request; check the published version every few seconds.

**D2. `/loci` stays on disk** (5.3 GB does not fit). It is already clustered by gene name, so one gene is a short sequential read. Multi-name searches merge per-name streams instead of sorting, and the concurrency limit applies.

**D3. Fixed pools and one shared annotation list.** A pool of 8 read-only connections per worker with 2 MiB caches replaces one connection per thread. The annotation list (about 7 MB) is loaded once per worker per version.

**D4. Responses built from pre-written pieces.** Each annotation's JSON is written once when the version loads. `/hits` joins the pieces for the matching annotations. Measured: 0.6 ms CPU for 2,500 matches, against about 31 ms for the current path, and a small memory peak (section 3).

**D5. Small tables that remove row-by-row counting.** A copy-count byte per annotation inside each posting, and a per-name total of symbol loci. `/names` then no longer scans `loci` rows, and `/hits` can report copies without touching `loci`.

**D6. Every endpoint is bounded.** A concurrency limit and a time budget on `/annotations/*` (today it has neither), and a polite `503` with `Retry-After` when busy. A summary table per annotation (namespace counts) so `/annotations/{id}/namespaces` does not open a multi-GB shard.

**D7. No spilling sorts, and temp files on NFS.** `/loci` multi-name terms use a per-name merge instead of `ORDER BY` over an `IN` list.

**D8. A real health and readiness check.** `/health` reports the published data version, that `presence.sqlite` and `loci.sqlite` open, that their versions match, and that a test lookup works. Today it always says `ok`.

**D9. Memory caps on child processes and cache headers.** `gffread` children get a data-size limit (128 MB) and are the first to be killed if memory is ever short. Responses carry an ETag and `Cache-Control` from the data version.

## 7. Work list

Order matters: D1 first, then the limits that make the budget true.

- [x] D1 index in RAM (`/dev/shm` copy, version check every few seconds)
- [x] Compose: `mem_limit`, `memswap_limit`, `shm_size`, log rotation, `TMPDIR` and `SQLITE_TMPDIR` on NFS
- [x] Caps: `HITS_MAX_INFLIGHT=3`, `LOCI_MAX_INFLIGHT=3`, `SEQUENCES_MAX_INFLIGHT=1`, SQLite caches, `--max-requests`, `SEQUENCES_MAX_INDEX_MB=256`
- [x] D3 connection cap (thread pool of 12 per worker) and shared annotation list
- [x] D6 limits and time budgets on `/annotations/*`
- [x] D8 real health check
- [x] `tools/check_api_memory.py`: replay worst-case requests, assert peak memory per the table in section 4
- [ ] D4 pre-built responses (then raise `HITS_MAX_INFLIGHT`)
- [ ] D5 count tables, D7 per-name merge, D9 child limits and ETag

## 8. Progress log

- 2026-11-24: decisions locked; worst-case `/hits` memory measured (section 3); budget drafted.
- 2026-11-24: D1 done in code (`helpers/presence.py`: copy to `PRESENCE_RAM_DIR` under a file lock, fallback to disk when too large / no space, no network-drive calls between checks every `PRESENCE_RECHECK_S`), 5 new tests in `tests/test_presence_ram.py`. Compose now has `mem_limit`/`memswap_limit` 2g, `shm_size` 1g, log rotation, scratch dir on NFS. Test run: 349 pass; `test_committed_openapi_schema_matches_the_app` fails with or without this change (checked by stashing it).
- 2026-11-24: caps done as defaults in `api/settings.py` (HITS 3, LOCI 3, SEQUENCES 1, SQLite cache 2 MiB, `SEQUENCES_MAX_INDEX_MB` 256), SQLite `soft_heap_limit` 64 MiB (presence and legacy corpus reader), gunicorn `--max-requests 5000 --max-requests-jitter 500`. Tests: 349 pass (openapi test excluded, known failure). Still open from the settings table: gate on `/annotations/*` (D6) and the connection pool (D3).
- 2026-11-24: D3 and D6 done. `/annotations/*` has a gate (3 per worker, 503 + `Retry-After: 1` when full) and a 20 s query budget per shard (`ANNOTATIONS_QUERY_BUDGET_S`, 503 `query_too_slow`). Handler threads cut from 40 to 12 per worker, which caps open connections; the annotation list is now one copy per version instead of one per thread (was up to 40 x 7 MB per worker). Budget table adjusted (connections 24 MB per worker instead of 16). Tests: 353 pass (openapi test excluded).
- 2026-11-24: D8 done: `/health` now opens the published index, runs a test lookup, checks `loci.sqlite` against the index version, and reports `data_version`, `index` (ram or disk), `lookup`, `loci`; returns 503 `unavailable` with a reason when nothing is published or versions differ. Without `PRESENCE_DIR` it still answers plain `ok`. `tools/check_api_memory.py` written and run: 144 worst-case requests (50 names present in 18,000 annotations, compact, records and TSV) from 24 client threads against the real app: 51 answered, 93 refused with 503 by the gates; peak RSS 203 MB for one worker (client side included) against the 486 MB limit. `/loci` is not part of this replay (needs a loci index); its memory is bounded by the 3-slot gate and 16 MiB cache. Tests: 355 pass (openapi test excluded).
