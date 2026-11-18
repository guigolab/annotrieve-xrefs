# Build plan: per-annotation sidecars, manifest and the presence table

Status: **built and run on the real corpus** (sidecars, manifest, `presence.sqlite`, `loci.sqlite`, `loci_ref.sqlite`; results in section 8). Sync with new/removed annotations (D6) and namespaces (D7) are still open. Read `API-DESIGN.md` (endpoints) and `PLAN-serve-index-and-cleanup.md` (served files, memory, cleanup) first. This document is only about **the script and logic that build the new tables**.

## 1. Scope

Decided:

- **D1:** gene coordinates go into the sidecar now (one pass over each shard serves both `presence.sqlite` and the later `loci.sqlite`).
- **D3:** a **new manifest** is the source of annotation ordinals and metadata.
- **D5:** shards are **rewritten in place**, so sidecars must detect staleness by fingerprint, not by path.
- **Sources:** an allow-list: only `source_database` Ensembl, GenBank and RefSeq are indexed. Community annotations are already skipped upstream by the import; if one shows up anyway (a `source_database` containing `community`) it is recorded `skipped_community`, and any other unknown value `failed`. Neither gets a sidecar.
- **Namespaces (Pfam, GO, ...) are deferred entirely.** No namespace aggregate, no `namespaces.json`, no Annotrieve filter work in this plan. The sidecar format reserves room for it (section 3.4), nothing more.
- **D6 (sync of new/replaced/removed annotations) and D7: postponed.**
- **No validation against the old global tables.** The corpus has a different number of annotations than when `xref_file`, `gene_hit` and `gene_xref` were built, so they are not a reference any more. Validation is against the shards themselves (section 7).

In scope now:

1. A **per-annotation sidecar** (names and gene coordinates) written from each shard.
2. A **manifest** (`manifest.sqlite`): one row per annotation, ordinals, metadata and fingerprints.
3. **`presence.sqlite` built from the sidecars and the manifest.**

Later: `loci.sqlite` / `/loci` (its input already exists in the sidecar), the sync mechanism, namespaces.

## 2. What the shard gives us

Per-annotation shard `per_gff/<annotation_id>/genes.sqlite` (`gene_corpus/sql/per_gff.py`):

| Table | Used for |
|---|---|
| `gene` | `local_id` (PK), `source_gene_id`, `feature_type`, `seqid`, `start`, `end`, `strand`, `biotype`, `primary_name` |
| `gene_xref` | rows with `namespace` `symbol` or `alias`: `local_id`, `accession` (the name). Indexes `xref_resolve(namespace, accession, local_id)` and `xref_by_gene(local_id, namespace)` |
| `meta` | `annotation_id`, `taxid`, `assembly_accession`, `organism_name`, `source_database`, `source_provider`, `profile_id`, `gff_path` |

Measured earlier (Phase 0, an older corpus snapshot, so order of magnitude only): about 18k shards, 418.7 GiB in total, median 14.7 MiB, max 3.5 GiB; NFS cold 16 KiB random read p50 19 ms; sequential about 58 MiB/s.

Not in the shard `meta`: `assembly_name`, `release_date`, and the taxon lineage. See D8.

## 3. Artifact A: the per-annotation sidecar

### 3.1 Purpose
One small file per annotation with everything the presence build needs, and everything the later loci build needs, so the weekly build never opens the shards or the global DB.

### 3.2 Location and lifecycle
- `per_gff/<annotation_id>/sidecar.tsv.gz`, next to `genes.sqlite` (decided); atomic write (`.tmp`, then rename). Gzip from the standard library.
- **Fingerprint:** because shards are rewritten in place (D5), the header stores the shard's `size`, `mtime_ns` and the `meta` values, plus a hash of the shard's `gene` row count and `max(local_id)`. A sidecar is stale when any of them differs from the current shard. Cheap to check (a `stat` and two tiny queries).
- Removing an annotation directory removes its sidecar.

### 3.3 Content (sidecar_version 1)
Header lines (`#key=value`): `sidecar_version`, `annotation_id`, `n_genes`, `n_name_rows`, the fingerprint fields, `built_at`.

Section `#names`: one line per distinct `(name, kind)`, sorted by `name` (byte order), then `kind`:

```
name<TAB>kind<TAB>local_id
```

- `name`: `gene_xref.accession`, **as stored in the shard**. Whether it is already casefolded is D4; the writer casefolds it (`str.casefold()`) either way, and records in the header whether any value changed.
- `kind`: `s` (`symbol`) or `a` (`alias`).
- One line per `(name, kind, local_id)`: a name that belongs to two genes of the same annotation yields two lines, which is what the loci build needs; presence collapses them.

Section `#genes`: one line per **named** gene (a gene with at least one `#names` row; decision D10), sorted by `local_id`. Unnamed genes are not written: they cannot be reached by name and, on the giant unannotated shards, dominated time and size. The header keeps `n_genes` (all genes in the shard) and adds `n_genes_named`:

```
local_id<TAB>source_gene_id<TAB>seqid<TAB>start<TAB>end<TAB>strand<TAB>feature_type<TAB>biotype<TAB>primary_name
```

Because `#names` references `local_id` only, coordinates are stored once per gene, not once per name.

### 3.4 Reserved for later
A `#namespaces` section may be added in a later `sidecar_version` once the namespace logic is decided. Nothing is written for it now.

### 3.5 Size
Unmeasured. The earlier guess was 100-150 KB for names only; adding genes makes it larger. M2 measures it.

## 4. Artifact B: the manifest

### 4.1 Purpose
The single record of what is indexed and under which ordinal. It replaces the dependence on `annotation_key` in the global DB (D3) and is the hook the later sync (D6) will use.

### 4.2 Schema (`manifest.sqlite`, small)

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);   -- manifest_version, created_at

CREATE TABLE annotation (
    ord INTEGER PRIMARY KEY,              -- assigned once, never reused
    annotation_id TEXT NOT NULL UNIQUE,
    taxid INTEGER,
    organism_name TEXT,
    assembly_accession TEXT,
    assembly_name TEXT,                   -- source: D8
    source_database TEXT NOT NULL,        -- Ensembl | GenBank | RefSeq
    source_provider TEXT,
    release_date TEXT,                    -- source: D8
    gff_path TEXT,
    shard_fingerprint TEXT,               -- as seen when the sidecar was written
    sidecar_fingerprint TEXT,             -- hash of the sidecar bytes
    n_genes INTEGER,
    n_name_rows INTEGER,
    status TEXT NOT NULL                  -- indexed | pending | failed | skipped_community | removed
);
CREATE INDEX annotation_status ON annotation(status);

CREATE TABLE lineage (taxid INTEGER NOT NULL, ord INTEGER NOT NULL, PRIMARY KEY (taxid, ord)) WITHOUT ROWID;  -- source: D8
```

### 4.3 Population (for now)
`serve_build sidecars` walks `per_gff/*`, reads each shard's `meta`, and upserts a manifest row (new `ord` only when `annotation_id` is new). It does **not** reconcile against Annotrieve: that is D6. A shard that is no longer on disk is left as it is in the manifest until D6 defines removal.

## 5. Artifact C: `presence.sqlite` from sidecars

Output schema is unchanged from `PLAN-serve-index-and-cleanup.md` section 3.1 (`meta`, `annotation`, `lineage`, `name_dict`, `posting`); the posting blob codec is the existing `gene_corpus/sql/postings.py`. Inputs are the manifest and the sidecars with status `indexed`.

### 5.1 Algorithm
1. Read `annotation` and `lineage` from the manifest; keep only `indexed` rows from allowed sources.
2. Level 1: for batches of about 256 sidecars (file-descriptor limit), read only the `#names` section, drop `local_id`, and write sorted intermediate runs of `(name, kind, ord)`. Runs go to `--tmp-dir`, which must be on the NFS work directory, **not local disk**.
3. Level 2: k-way merge of the runs (`heapq.merge`). Rows with the same name form a group; within a group each kind's ords are strictly increasing.
4. For each name in sorted order: assign `name_id`, compute `n_sym` and `n_alias_only`, encode `sym` and `alias` blobs, write `name_dict` and `posting`.
5. Write `meta` (`data_version`, `schema_version`, counts, hash of the manifest state) and publish through the existing versioned-directory and atomic `current` symlink (`publish_version`).

Memory is bounded by the batch size and one name group (a very common name has about 18k ords).

### 5.2 Reuse
`postings.py`, the `presence.sqlite` writer and its checks, `publish_version`, and the tests are reused. The existing `xref_file` input path is **legacy**: it is kept only until the sidecar path replaces it, and is not used as a reference (section 7).

## 6. CLI

Implemented (M3), all under `python -m gene_corpus.serve_build`:

```
catalogue --work-dir DIR [--api-base URL | --url URL | --tsv-file F]   # report TSV -> manifest (D8)
sidecars  --work-dir DIR [--workers 2] [--force] [--ids FILE] [--sample N] [--limit N]
status    --work-dir DIR
```

- Defaults: shards in `<work-dir>/per_gff`, manifest at `<work-dir>/serve/manifest.sqlite`, sidecars next to the shards.
- `catalogue` streams the report (flat memory), replaces the manifest `catalogue` table, fills `assembly_name`, `release_date`, `in_catalogue` and the `lineage` table (each annotation's taxon plus its `taxon_lineage`), and prints shard-only / catalogue-only / database-mismatch counts. It removes nothing (D6). A bad report keeps the previous catalogue.
- `sidecars` runs shards in parallel but only the main process writes the manifest. Results are consumed in sorted shard order, so ordinals are deterministic. It skips fresh sidecars, rewrites stale ones, records `skipped_community` and `failed` shards without stopping, commits every 200 shards, saves progress on Ctrl-C, prints rate and ETA, and re-applies a stored catalogue at the end.
- `status` prints counts per status, catalogue coverage and the number of indexed annotations with no names.
- `build` (M4): `serve_build build --work-dir DIR [--name-prefix tp5] [--limit-annotations N] [--batch 256]` merges the sidecars into `presence.sqlite` (two-level merge, runs in `<work-dir>/serve/tmp`). Includes only `indexed` and `in_catalogue=1` annotations. Partial builds go to `partial-<version>/` and never switch `current`.
- Also: `loci` (`loci.sqlite`) and `loci-ref` (`loci_ref.sqlite`), see section 8. Not implemented: a separate `verify` command; `build` and `loci` verify as they go. The old `build` from `xref_file` (no sub-command) is legacy.

## 7. Validation (no old global tables)

1. **Sidecar vs shard:** for sampled shards, the sidecar `#names` set equals a fresh `SELECT` of `symbol`/`alias` rows from `gene_xref` (casefolded), and `#genes` equals `SELECT ... FROM gene`; `n_genes` equals `COUNT(*)`.
2. **Fingerprint:** modifying a shard in place makes `sidecars` rewrite exactly that sidecar.
3. **Presence vs sidecars (self-consistency):** sum of `n_sym` equals the number of distinct `(name, ord)` symbol pairs counted by an independent pass over the sidecars; for sampled names, the decoded postings equal a direct scan of the sidecars; every manifest `indexed` annotation appears in at least one posting (unless it has no names, which is reported).
4. **Brute-force query check:** for about 2,000 random terms (exact and prefix), the annotation set computed from `name_dict` and `posting` equals the set computed by scanning the sidecars directly.
5. **Build report:** counts at each step, skipped (`Community`), failed, and no-name annotations, timings.
6. Anything that looks odd is reported, not fixed silently. Corpus-level questions about old tables (such as the 21,563 pairs seen only in `xref_file`) are dropped from this plan because those tables are no longer a reference.

## 8. Milestones and gates

| M | Work | Exit / gate |
|---|---|---|
| M1 | Sidecar writer and reader for one shard, format spec, manifest schema, unit tests (round trip, determinism, atomic write, fingerprint and stale detection, D4 casefold flag) | tests pass |
| M2 | **Measure on about 200 random shards** (cold NFS): time per shard, sidecar bytes (names plus genes), projected full backfill time and total size. Also report source_database/provider values found, and the D4 answer | report; **stop and re-decide if the projection is bad** |
| M3 | `catalogue` command: stream the report TSV (D8) into the manifest (`assembly_name`, `release_date`, `lineage`). Parallel backfill `sidecars` (`--workers`, resumable, sidecars next to the shards at `per_gff/<id>/sidecar.tsv.gz`). Report shards without a catalogue row and catalogue rows without a shard (listed only, no removal: D6) | every allowed shard has a sidecar; manifest metadata filled; mismatch lists reviewed |
| M4 (done) | Merge builder writing `presence.sqlite`; trial on a name prefix (for example `tp5`) | validation steps 3-4 pass on the trial |
| M5 (done 2026-11-18) | Full build, report and timings | validation clean or every difference explained |
| M6 (done: `PRESENCE_DIR`, alias as suggestion, live; a formal query replay was not run) | Publish, API reader switch behind a flag, query replay | zero unexplained diffs |

### M2 results (205 shards: 5 largest + 200 random, 2 workers, cold-ish NFS)

First run (all genes in `#genes`) and second run (D10: named genes only), same sample:

| Fact | All genes | Named genes only (current) |
|---|---|---|
| Peak writer memory | 51 MiB | **51 MiB** |
| Validation vs fresh shard query | 0 failures | **0 failures of 205** |
| Wall time for the 205 shards, 2 workers | 10.0 min | **59 s** |
| Time per shard mean / p50 / p95 / max | 4.0 s / 0.6 s / 1.7 s* / 3.3 min | **0.5 s / 0.3 s / 1.7 s / 2.4 s** |
| Projected single-process backfill (all 17,984) | 3.4-5.7 h | **39 min (by bytes) to 2.3 h (by mean)** |
| Sidecar / shard size ratio | 0.033 | **0.0022** |
| Projected total sidecar size | about 14 GiB | **2.3 GiB measured in the full backfill** (the earlier 0.9 GiB projection was wrong) |
| Name rows in the sample | 1,471,877 | 1,471,877 (identical: nothing lost) |

\* first-run p95 was 2.4 s.

Other facts: names are already casefolded in all 205 shards (D4); `source_database` was only GenBank, Ensembl and RefSeq (D9). The giant shards (4.6-9.4M genes) have no symbol/alias rows, so their sidecars are about 400 bytes and take about 2 s (they were 100-200 s). In Ensembl shards about half the genes are named (for example 43k of 80k). Gate **passed**.

### M3 results (full backfill, 2026-11-17)

17,984 shards indexed, 0 failed, 84.4 min with 3 workers (3.6 shards/s), peak worker RSS 70 MiB, 2,361 MiB of sidecars. 16,578 shards are in the live catalogue, 1,406 are not (excluded from presence); the catalogue has 19,238 rows, so 2,660 annotations have no shard yet (D6). 7,837 indexed shards (43.6%) have no searchable names.

### M3b: shards for catalogue annotations that had none (2026-11-18)

`serve_build missing` found 2,660 catalogue annotations without a shard: 2,039 CommunityRegistry (skipped), 621 harvestable (GenBank 506, Ensembl 69, RefSeq 46). Harvested with `harvest_cli --annotation-ids-file missing_ids.txt --files-root /data/annotrieve/files --workers 2`: 619 ok in the first pass (peak 89 MiB for 5-annotation pilot); 2 large annotations failed under `ulimit -v 1500000` (empty `MemoryError` message) and succeeded with `--memory-flush-units 200000` (23 min, peak 554 MiB). For big annotations use a small flush window. Then `sidecars --ids`, `catalogue`, `build`.

### M5 results (first full presence build, 2026-11-18)

`serve_build build`: 17,199 annotations (`indexed` and `in_catalogue`), 9,545 of them with names; 3,568,990 distinct names (1,524,416 alias-only); 66.1M symbol pairs, 2.28M alias pairs; 38 runs; `presence.sqlite` 310 MB (postings 72 MB; 3.65M varint, 10,354 bitset blobs). 31 min (level 2 about 25 min), **peak RSS 334 MiB**. Verification: 31 annotations, 233,073 name rows, all found. Published as `serve/v9e8dc2f914f42ba5`, `current` points at it.

Smoke test through the API service code (`tools/smoke_presence.py`, cold NFS): open reader 116 ms; `/hits tp53` 81 ms (2,492 symbol annotations, 2 alias-only extra); `/hits` for 3 names 74 ms; Mammalia scope 9 ms; `/names` prefix 0.1-7.5 ms; process peak 30 MiB. Known issue: most alias-only names are locus-style IDs (`a1f95_gene3136`, `afu2g05510`) and some contain a comma (`agap_ensangg...,ensangg...`), which cannot be sent in a comma-separated `names=` list.

### Loci build (done; measured 2026-11-18)

`serve_build loci --work-dir DIR` writes `loci.sqlite` next to the published `presence.sqlite` (same `v<data_version>/`), from the sidecars' `#names` joined to `#genes` per annotation (batches by named-gene budget, 300,000 genes in memory, larger annotations use a scratch SQLite file under `serve/tmp`; two-level merge as in the presence build). `locus(name_id, ord, local_id, seqid_id, start, end, strand, feature_type_id, biotype_id, via, gene_id, primary_name)` clustered by `name_id`, plus `seqid`, dictionaries and `name_loci(name_id, n_loci)`. Decision: rows carry the source gene id; `locus_ref = annotation_id:local_id`. The build walks `presence.name_dict` in lockstep and fails on any difference in the annotation set, the name list or the symbol/alias counts per name, so the two files cannot disagree. The real run produced the file now served by `/loci` (80,995,350 genes in `loci_ref`).

`serve_build loci-ref --work-dir DIR` (also run at the end of `loci`) writes `loci_ref.sqlite` next to `loci.sqlite`: one row `(ord, local_id) -> name_id` per gene, so a `locus_ref` resolves to its `locus` row with the primary key. It is a separate file so the published `loci.sqlite` is never modified; it is built by one sorted `GROUP BY` over `loci.sqlite` (sort spills to `serve/tmp` through `SQLITE_TMPDIR`, page cache 64 MiB). Run it alone on an already published version to add the file without rebuilding `loci.sqlite`. Measured on the real index (2026-11-18): 80,995,350 genes, 1.31 GB, 134 s, peak RSS 157 MiB.

## 9. Open decisions

| # | Question | Needed by |
|---|---|---|
| D4 | **Answered by M2 (205 shards):** shard `symbol`/`alias` names are already casefolded (0 rows changed). The writer still casefolds as a safety net. | done |
| D8 | **Decided:** metadata comes from the Annotrieve report endpoint `GET /annotations/report?selected_fields=release_date,taxon_lineage` (TSV, streamed line by line). Columns: `annotation_id, assembly_accession, assembly_name, organism_name, taxid, database, provider, source_url, bgzip_path, csi_path, release_date, taxon_lineage` (`taxon_lineage` is `;`-joined taxids; dates ISO). Joined to shards by `annotation_id`. Base URL and any API-side cap are checked in M3 | M3 |
| D9 | **Decided:** the manifest uses an allow-list (`Ensembl`, `GenBank`, `RefSeq`, casefolded); any `source_database` containing `community` is recorded `skipped_community`, any other unknown value `failed`; neither gets a sidecar. Community GFFs are already skipped upstream by the import, which is why M2 saw none. Facts from M2: in 205 shards `source_database` is only `GenBank` (106), `Ensembl` (74), `RefSeq` (25); **no `Community` value appeared**. `source_provider` has 16 shards with `community` (lowercase) and many institutions. Which field marks "Community"? Open | before M3 |
| D10 | **Decided:** `#genes` holds only genes with at least one symbol/alias row. Re-measure with M2 (same seed) to see the new time and size | before M3 |
| D6 | Sync with new, replaced and removed annotations | postponed |
| D7 | Namespaces and the Annotrieve annotation filter | postponed |
| D2 | Which namespaces | postponed with D7 |

## 10. Risks

- **Harvest name gate: strict, kept on purpose (decision 2026-11-18).** `stream/symbol_gate.py` drops names matching `[A-Za-z]+_\d+$` (and `LOC\d+`, `_LOCUS`, `gene-123`, and echoes of `locus_tag`/`ID`). This also drops real copy-suffixed symbols (`Polq_1`, `Lrba_4`, `selenot_1`; `selenot2_1` passes because of its digit). Measured on 40 named shards: 1,665 of 519,242 genes without a symbol (0.3%). A relaxed rule was tried and measured over all 16,578 shards (`tools/measure_gate_relax.py`): 526 shards (3.2%) and 2,581,872 genes would gain names, but 91% of them in 10 shards and mostly locus or transposon IDs (`TE_00001670`, `LMJFC_010005100`, `SRAE_1000000100`, `EAI_06374`); about 500 other shards would gain real copy-suffixed symbols. The risk of flooding the index with locus genes outweighs that loss, so the gate was reverted and no shard is re-harvested. If revisited: reject letters + `_` + 4 or more digits (`\d{4,}`) for symbols only, aliases stay strict, then re-harvest the shards in `gate_gain_full.tsv` with `harvest_cli --annotation-ids-file`. The 7,837 nameless shards are genuinely nameless (checked on 40 samples; locus tags only). Locus tags are not searchable in v1 (`locus_tag` namespace, D2/D7). Whether old shards were harvested with a looser gate is unchecked. Scripts: `tools/diagnose_nameless.py`, `tools/measure_gate_relax.py`.

- **Backfill cost over NFS is unknown.** M2 decides. If too slow, we cut what the sidecar contains or parallelize harder; we do not go back to the old global tables.
- **Sidecar size is unknown** with coordinates included; the weekly merge reads only `#names`, so it does not pay for `#genes`, but disk and the backfill do.
- **In-place shard rewrites** (D5) mean a fingerprint bug could serve stale data. The fingerprint test in M1 is not optional.
- **Many small files on NFS:** 18k sidecars plus writes; batch per worker.
