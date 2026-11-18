# Serve files, build pipeline and cleanup

Status: **served files built and live** (`presence.sqlite`, `loci.sqlite`, `loci_ref.sqlite`; endpoints in `API-DESIGN.md`). Still open: sync (D6), the cleanup in section 6, namespaces (D7). **How the tables are built is in `PLAN-build-tables.md`**; this document keeps the served-file definitions, serving and memory, and the cleanup lists. Earlier versions assumed DFS-ordered ordinals, a delta file and a local-disk copy. All three are dropped; see "What changed" below.

Remaining: the sync mechanism for new, replaced and removed annotations, and the cleanup in section 6.

## 1. What was measured (Phase 0, real corpus)

| Fact | Value |
|---|---|
| Annotations | 17,984 (Ensembl 5,988, GenBank 9,428, RefSeq 2,568), keys contiguous 1-17,984 |
| `gene_corpus.sqlite` | **82.1 GiB** on NFS 4.2 (rsize 1 MiB, `hard`) |
| Rows in the big tables (`sqlite_stat1`) | `gene_hit` 559M, `xref_meta` 244M; `xref_file` includes about 80M symbol/alias rows |
| Symbol rows / alias rows | 77.7M / 2.25M |
| Distinct symbol/alias names | 3.93M (71% occur in one annotation only; 18.6k names occur in more than 1,000) |
| Lineage | 360k rows, 19,349 taxids |
| Per-annotation shards | 17,984 files, **418.7 GiB** (median 14.7 MiB, max 3.5 GiB) |
| NFS cold reads | 16 KiB random read p50 19 ms, p95 37 ms; sequential 58 MiB/s |
| Current SQL, warm | 10 very common genes about 0.1 s; 200 very common genes (1.25M rows) 1.9-2.6 s; broad `/names` prefixes 0.3-0.6 s |
| Live API, cold | 10 genes 2.8 s; `/names?prefix=sel` 4.9 s; `hox` 7.6 s (ends in `prefix_too_broad`) |
| Estimated compact name index | about 183 MiB (about 230 MiB with page slack) |

Reading: the slow cases are **cold random reads on a file far larger than RAM**. Warm SQL is fine. The fix is to serve from small, name-clustered files that the page cache can hold or read sequentially.

### Phase 0b results (gene_hit and consistency)

- **Consistency:** 30 random windows per namespace (about 1.15M pairs each for symbol and alias) show **0 mismatches** between `xref_file` and `gene_hit`. A dry run found 21,563 pairs only in `xref_file` in the very first alphabetical slice; not seen elsewhere (inspect it separately).
- **Reading loci from `gene_hit` by name:** common names (3k-13k rows) took 0.07-0.25 s in the second run, 0.02-0.04 s when repeated, and up to about 2.5 s the first time they were ever read (dry run). Small names took 1-13 ms. Cold cost depends on NFS cache state, so the true worst case is about 2.5 s per common name; typical is far lower.
- **Average row payload** in `gene_hit` is about 54 bytes (text fields plus integers).
- **Case:** 77% of symbol rows have `primary_name` equal to the accession ignoring case, 18% identical, 4.6% a different name (such as an Ensembl id). 86% of alias rows have a different primary name.

## 2. What changed from the earlier plan

- **Dropped DFS ordinals.** `annotation_key` already grows monotonically, so use it directly as the ordinal. Taxon scope comes from the small `annotation_lineage` table.
- **Dropped the delta file.** The derived files are swapped in atomically, so the API stays online. The earlier estimate "loci file about 10-15 min" was wrong: streaming `gene_hit` cold takes hours. The build source is now **per-annotation sidecars** (see `PLAN-build-tables.md`); the old `xref_file` build is legacy and **not a reference**: the corpus has changed since those tables were built, so validation is against the shards and sidecars.
- **Dropped the local-disk copy.** You do not want to use local disk. The presence file is small enough to live in the kernel page cache after one warm-up read.
- **Added a second file for coordinates** (`loci`), because `/loci` is now a core endpoint.

## 3. The two served files

Both are published as a **pair** with the same `data_version` once both exist; until `loci.sqlite` is built, `presence.sqlite` is published alone. The API refuses to start with mismatched versions of files that are present. Build inputs: per-annotation sidecars (`PLAN-build-tables.md`), plus the small `annotation` and `annotation_lineage` tables.

### 3.1 `presence.sqlite` (about 230-300 MiB): endpoints `/names`, `/hits`

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
-- data_version, schema_version, built_at, n_annotations, source_db_mtime

CREATE TABLE annotation (            -- about 18k rows, loaded into memory at startup
    ord INTEGER PRIMARY KEY,         -- = annotation_key
    annotation_id TEXT NOT NULL UNIQUE,
    taxid INTEGER NOT NULL,
    organism_name TEXT,
    assembly_accession TEXT,
    assembly_name TEXT,
    source INTEGER NOT NULL,         -- 1 Ensembl, 2 GenBank, 3 RefSeq (helpers/source_codes.py)
    release_date TEXT,
    gff_path TEXT                    -- for download_url / gff links
);

CREATE TABLE lineage (               -- copy of annotation_lineage (360k rows)
    taxid INTEGER NOT NULL,
    ord INTEGER NOT NULL,
    PRIMARY KEY (taxid, ord)
) WITHOUT ROWID;                      -- scope = one range read per taxid, cached in an LRU

CREATE TABLE name_dict (             -- 3.9M rows; THE table that prefix expansion scans
    name TEXT PRIMARY KEY,           -- casefolded symbol or alias
    name_id INTEGER NOT NULL UNIQUE,
    n_sym INTEGER NOT NULL,          -- distinct annotations where the name is a symbol
    n_alias_only INTEGER NOT NULL,   -- annotations where it is only an alias (the 'suggestion' count)
    -- (n_loci and category are NOT stored in v1: n_loci needs gene_hit and arrives with the loci build;
    --  the API classifies the at most 500 expanded names at request time, since the classifier lives in api/)
) WITHOUT ROWID;

CREATE TABLE posting (               -- blobs, read only by name_id
    name_id INTEGER PRIMARY KEY,
    sym BLOB,                        -- namespace=symbol postings (NULL if none)
    alias BLOB                       -- namespace=alias postings (NULL if none)
);
```

**Alias semantics:** `sym` and `alias` are separate blobs. A match is `sym`; the suggestion for a term is `alias AND NOT sym`. Names that exist only as aliases (`n_sym = 0`) are hidden from prefix expansion unless `include_alias_only=true`.

Why split `name_dict` and `posting`: **prefix expansion** scans `name_dict` by range (up to 501 rows). Keeping blobs out of that table means the scan reads dense small rows (about 100 MB in total, always cached), and the blobs are read only for the names actually used. `n_sym` and `n_alias_only` make `prefix_too_broad` refusals free; `n_loci` will be added with the loci build.

Posting blob: first byte is a format tag. `0x01` = sorted ordinals as LEB128 varint deltas. `0x02` = raw bitset of `ceil(N/8)` bytes (used for dense names, about 2.3 KB at N=17,984). Pick the smaller. Python ints serve as bitsets: prefix unions are `OR` over the names' bitsets, which gives exact union counts.

Query flow for `/hits`:
1. Expand each term (exact: one `name_dict` lookup; prefix: range scan, minus `exclude`, filtered by `categories`).
2. For each term, `OR` the bitsets of its names (symbol and alias parts kept separately for `matched_via`).
3. Scope mask: for each include taxid, read its `lineage` range (LRU cached); exclude subtracts. Apply `sources`.
4. Per annotation, collect matched term indices; apply `min_matched`; join to the in-memory `annotation` list.

### 3.2 `loci.sqlite` (guess 2.5-3.5 GB, to measure; **later**, after the presence table; its input, gene coordinates, is already in the sidecar): endpoint `/loci`, expansion for `/sequences`

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE locus (                 -- one row per (name, gene); clustered BY NAME
    name_id INTEGER NOT NULL,        -- same ids as presence.name_dict
    ord INTEGER NOT NULL,            -- annotation ordinal
    local_id INTEGER NOT NULL,       -- shard gene id (see API-DESIGN open decision 1)
    seqid_id INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand INTEGER NOT NULL,
    feature_type_id INTEGER NOT NULL,
    biotype_id INTEGER,
    via INTEGER NOT NULL,            -- 0 symbol, 1 alias
    primary_id INTEGER NOT NULL,     -- dict_display id of gene_hit.primary_name (the gene's own name, original case)
    PRIMARY KEY (name_id, ord, local_id)
) WITHOUT ROWID;

CREATE TABLE seqid (ord INTEGER NOT NULL, seqid_id INTEGER NOT NULL, seqid TEXT NOT NULL,
                    PRIMARY KEY (ord, seqid_id)) WITHOUT ROWID;
CREATE TABLE dict_feature_type (id INTEGER PRIMARY KEY, name TEXT UNIQUE);
CREATE TABLE dict_biotype (id INTEGER PRIMARY KEY, name TEXT UNIQUE);
CREATE TABLE dict_display (id INTEGER PRIMARY KEY, name TEXT UNIQUE);  -- original-case primary names (TP53, Trp53, Oprl1); from gene_hit.primary_name
-- Display rule: if casefold(primary name) == matched name, show the primary name; otherwise show the matched name (lowercase).
-- Phase 0b: 77% of symbol rows differ from the accession only in case, 4.6% have a different primary name (e.g. an Ensembl id), and 86% of alias rows have a different primary name.
```

Why clustered by name: all loci of one name sit in adjacent pages. A cold request costs a short sequential read (a very common gene is about 8k rows, a few hundred KB), not thousands of 19 ms random reads. A gene looked up through a separate `gene` table by `(ord, local_id)` would be random reads again, so coordinates are stored in the row.

Row size estimate: about 25-30 bytes, times 80M rows, is 2-3 GB. A guess. `gene_hit` has 559M rows for all namespaces, so this is a large reduction.

Query flow for `/loci`: expand terms (as above) to `name_id`s, range-scan `locus` per name, keep rows whose `ord` is in the scope mask, join `seqid` and `annotation` in batches, stream.

## 4. Build

**Superseded in part by `PLAN-build-tables.md`.** The target build source is now per-annotation sidecars merged k-way (section 5 of that document). Sections 4.2 to 4.4 below describe the working `xref_file` path, which is **legacy** (kept only until the sidecar path replaces it; not used as a validation reference); 4.3 (loci from `gene_hit`) is no longer the plan.

**Status: the presence build from `xref_file` is written** (`gene_corpus/sql/postings.py`, `gene_corpus/sql/serve_presence.py`, `gene_corpus/serve_build_cli.py`, tests in `gene_corpus/tests/test_serve_presence.py`, 16 passing). Run: `python -m gene_corpus.serve_build --work-dir /data/annotrieve/gene_corpus [--name-prefix tp5 --out-dir /tmp/trial]`. It has not been run on the real corpus. A synthetic 921k-row run built in 1.3 s, including a 200-name re-check; extrapolated, the 80M-row build is about 2 minutes of CPU plus NFS read time. Decisions applied: casefolded index with a `case_sensitive` filter on loci (Phase 0b: 25.7% of sampled names have more than one spelling, mostly species conventions such as ZFC3H1/Zfc3h1/zfc3h1); alias hits are suggestions, not matches.

### 4.1 Inputs and outputs
- In: `gene_corpus.sqlite` (read-only, `immutable=1`), `--out-dir`, `--mem-mb`, `--tmp-dir` (any scratch space; not required to be local, but fast).
- Out: a versioned directory `serve/v<data_version>/` containing `presence.sqlite` and `loci.sqlite`, validated, then published by switching a `serve/current` symlink in one atomic step (renaming two files one after the other would leave a window with mismatched versions). The API resolves the symlink target on each reopen; old readers keep the old files until they close.
- Cross-check built in: every (name, annotation) in the postings must have at least one locus row and vice versa; mismatches between `xref_file` and `gene_hit` are counted and written to a build report (the repo has `repair_*` scripts, which suggests the two tables have drifted before).

### 4.2 Presence build (from `xref_file`)
1. Stream `xref_file` rows for `namespace IN ('symbol','alias')` in primary-key order: `(namespace, accession, annotation_key)`. Symbol first, then alias, each ordered by accession. Merge the two streams by accession (the Phase 0 script already does this in about 31-220 s).
2. For each name: assign `name_id` in sorted order, count `n_sym` and `n_alias_only`, encode `sym` and `alias` postings, write `name_dict` and `posting`.
3. Copy `annotation` (plus `organism_name`, `assembly_name`, `release_date`, `gff_path`) and `annotation_lineage`; annotations referenced by `xref_file` but missing from `annotation` are dropped, as the old `/hits` path ignored them.
4. Memory is constant: both streams are ordered and consumed once.

### 4.2b Case sensitivity
- `xref_file` and `gene_hit.accession` are casefolded (harvest calls `casefold()`); `xref_file` cannot supply original-case names. The only original-case form is `gene.primary_name` (via `gene_hit`; the shard `gene_xref.display` column holds description prose, not the symbol).
- Default: index casefolded names; apply `case_sensitive=true` to `/loci` rows only.
- Exact-case index (keyed by exact spelling plus a fold lookup) would require building from `gene_hit` or the shards, estimated at 2-7 hours cold. Decide after the spelling check (`tools/phase0b_check.py`, section 4).

### 4.3 Loci build (from `gene_hit`) - **not the plan; kept for reference**
Streaming `gene_hit` cold takes hours (about 3k pairs/s in the dry run). The loci build will read coordinates from the per-annotation sidecars instead, because the sidecar already carries gene coordinates (decision D1, settled). Original sketch:
1. Stream `gene_hit` rows for `namespace IN ('symbol','alias')` in primary-key order (clustered by name, then annotation, then gene).
2. Look up `name_id` through a streaming merge against the sorted `name_dict` (the same sort order), so no in-memory dictionary of 3.9M names is needed.
3. Dictionary-encode `seqid`, `feature_type`, `biotype` per annotation. Write `locus` in key order, `alias_primary` for alias rows. Add `n_loci` into `name_dict` in a second pass or in the same pass if both builds run together.
4. `gene_hit` has no source gene id. See API-DESIGN open decision 1.

### 4.3b Estimated runtime of the presence build (estimate, not measured)
| Step | Estimate |
|---|---|
| Stream about 80M symbol/alias rows from `xref_file` | 1-4 min warm, up to about 4 min cold (Phase 0 scan: 31.5 s warm, 220 s cold) |
| Group by name and encode postings in Python | 1-2 min |
| Write 3.9M `name_dict` and `posting` rows | 0.5-1 min |
| Copy `annotation` and `lineage` | seconds |
| Validation | 1-2 min |
| **Total** | **about 5-8 min warm, 10-20 min cold**; memory stays small (streamed); output about 250 MB |

The main uncertainty is NFS cache state. Building from `gene_hit` instead would take hours.

### 4.4 Validation before publish
- Row counts: `name_dict`, sum of `n_sym` vs the source scan, sum of `n_loci` vs `locus` rows.
- 200 random names: decoded postings equal `xref_file` presence; loci equal `gene_hit` rows.
- Parity replay: about 2,000 queries (random names, example sets, with and without scope, prefix terms) against the old and new paths, comparing annotation sets, matched indices and counts. Cut over only on zero unexplained diffs.
- Write `data_version` (hash of source mtime plus row counts) to both `meta` tables.

### 4.5 When to run
After each Annotrieve import: import → taxonomy export → `sync` (harvest, merge, index) → `serve_build` → publish. Imports are weekly and vary in size. **How new, replaced and removed annotations reach the sidecars and the served files is postponed** until the first table works (`PLAN-build-tables.md` D6). Removed or replaced annotations must disappear from the served files; confirm how `sync` handles that today (open question 3).

## 5. Serving and memory (2-3 GB budget, NFS)

| Item | Setting |
|---|---|
| Workers | 2 gunicorn workers; each loads `annotation` and the lineage LRU (a few MB) |
| SQLite page cache | 16-32 MiB per connection (down from 256 MiB) |
| Threads | cap the threadpool at about 8 per worker |
| `presence.sqlite` | about 230-300 MiB; read once sequentially at startup and after each publish (`cat file > /dev/null` at 58 MiB/s takes a few seconds). It stays in the kernel page cache |
| `loci.sqlite` | 2-3 GB, not fully cacheable; access is clustered by name, so it works from cold with short sequential reads, and hot names stay cached |
| mmap | leave `mmap_size=0` on NFS |
| FASTA | `SEQUENCES_MAX_INFLIGHT=2`, 32 MB slice cache, gffread children: budget about 500 MB |

Caching in front: nginx key on the canonical query string, `ETag` from `data_version`, and a warm-up script after each publish.

## 6. Tables to drop

"Drop" means: stop building and publishing; remove after the new files are validated and one release has run on the new files. **Archive the old DB file before deleting anything.** The dropped tables hold non-symbol namespaces (GeneID, GO, Pfam, InterPro, orthologs) that no endpoint reads today; the shards still contain them, but rebuilding from the shards is slow (418 GiB).

### 6.1 Tables in `gene_corpus.sqlite` (82 GiB)
| Table | Fate | Reason |
|---|---|---|
| `gene_hit` (559M rows) | keep until the loci file is validated, then drop | replaced by `loci.sqlite` for the served namespaces; the other namespaces are not served in v1 |
| `xref_meta` (244M rows) | drop | not read by the API; counts now live in `name_dict` |
| `xref_file` | drop after cutover | replaced by `name_dict` plus `posting`. It is the build input until then |
| `xref_file_status`, `merge_done`, `namespace_stats` | drop | build bookkeeping |
| `xref_label`, `xref_label_alias`, `xref_label_fts*` | drop from the served DB (export first if GO/Pfam search is wanted later) | label FTS only matters for non-symbol namespaces |
| `annotation`, `annotation_lineage`, `meta` | keep | small; copied into `presence.sqlite` |

Because `gene_corpus.sqlite` stays the build input in v1, these tables are not deleted immediately. The point is that **the API no longer opens this file**. If you later want a smaller build input, the sources of truth are the shards, and `serve_build` could read `tier_a_counts` and the `gene` tables from them instead (slower: the shards total 418 GiB).

### 6.2 Code that becomes dead after the build input changes
Only once nothing reads `xref_*` / `gene_hit` apart from `serve_build`:

| Path | Reason |
|---|---|
| `gene_corpus/prune_xref*`, `sql/prune_xref_file.py` | maintains `xref_file` |
| `gene_corpus/repair*`, `sql/repair_quotes.py`, `repair_gene_hit_deletes*`, `repair_symbol_deletes*` | repairs of the global tables |
| `gene_corpus/audit_*.py`, `selenoprotein_xref_counts.py` | ad-hoc audits on the global tables |
| `gene_corpus/labels/`, `namespaces/labels.py`, label tables in `sql/schema.py` | non-symbol namespaces (park, do not delete, if GO/Pfam search may return) |
| their tests | follow the code |

### 6.3 Keep
Harvest (`harvest_cli.py`, `stream/`, `profiles/`, `gating.py`, `resolve.py`, `catalog/`), `sql/per_gff.py`, `sql/lineage.py`, `seed_release_date`, `sync`, `merge`, `index` (they still feed `serve_build`), `namespaces/tiers.py`.

Do **not** remove `merge` / `index` / `attach` in v1: they are what keeps `xref_file` and `gene_hit` current. They go only if `serve_build` is later changed to read from the shards directly.

### 6.4 API code that changes
- `helpers/session.py`: open the two files, smaller caches; check that the `data_version` of both matches.
- `services/hits.py`, `services/names.py`: rewrite against `name_dict`, `posting` and the lineage LRU. Remove chunked `IN` helpers and the `annotation_lineage` joins.
- New `services/loci.py` and `routers/loci.py`.
- `services/sequences.py`: add `query` mode (expand via `locus`), `estimate`, a global NCBI limiter, then jobs.
- `helpers/names.py`: accept a trailing `*` as a prefix term.

## 7. Phases (staged: each stage ships on its own)

| Stage | Work | New file? | Exit |
|---|---|---|---|
| 1 API shape (done, on the new indexes instead of today's tables) | `/hits`: records, TSV, `min_matched`, `sources`, prefix terms (expansion by the existing `xref_file` range scan). `/loci` reads `gene_hit` by name (json, tsv, bed, `dry_run`, `gff_url`, `locus_ref = annotation_id:local_id`, display rule from section 3.2). `/sequences`: `query` mode, `estimate`, global NCBI limiter, `transcript=longest`, `on_error=skip`. Config: thread cap, 16-32 MiB per-connection cache, canonical nginx cache key, cheap `prefix_too_broad` from `xref_meta` counts | no | UI can be built against the final shapes; real latencies collected from logs |
| 2 Presence file (done) | Sidecars, manifest, then `presence.sqlite` from sidecars (`PLAN-build-tables.md` M1-M6), query replay against the sidecars, switch `/names` and `/hits` behind `HITS_BACKEND` | presence | Cold `/hits` and broad prefixes in milliseconds after warm-up; replay clean |
| 2b Sync (**next**) | Decide how new, replaced and removed annotations update sidecars and served files | no | Decision recorded (D6) |
| 3 Loci file (done) | Build `loci.sqlite` from sidecars, switch `/loci` | loci | Cold common-name read well under a second; build time acceptable |
| 4 FASTA jobs (done) | Async jobs, manifest, result cache | no | Large exports complete without a long-held request |
| 5 UI | Minor changes to the existing Distribution, Annotations and Compare views: one query object, loci TSV export, estimate and jobs in the FASTA download (see `API-DESIGN.md` section 7) | no | |
| 6 Cleanup | Section 6 | no | Smaller repo and published artifacts |

## 8. Open questions

1. `gene_hit` has no source gene id (API-DESIGN open decision 1). Which option?
2. How many symbol/alias rows does `gene_hit` hold? It should equal `xref_file`'s 80M if the two stayed in step; confirm.
3. How does `sync` remove deleted or replaced annotations from `xref_file` and `gene_hit`? If it does not, the build must filter against the live catalogue.
4. Annotation download URL pattern and the slice offset (API-DESIGN open decision 2).
5. Should `serve_build` stay a build-from-global-DB step, or eventually read the shards directly (and let `merge`/`index` retire)? The sidecar plan moves toward the second; the open decisions are in `PLAN-build-tables.md` section 9.
6. Namespaces (Pfam, GO, ...) and the Annotrieve annotation filter are deferred (`PLAN-build-tables.md` D7). Only `Ensembl`, `GenBank` and `RefSeq` annotations are indexed; `Community` is skipped.
