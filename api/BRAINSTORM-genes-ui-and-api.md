# Genes page and xrefs API: brainstorm

> **Status: historical brainstorm.** The agreed direction is now in `../API-DESIGN.md` (three endpoints: `/hits`, `/loci`, `/sequences`, with prefix terms) and `../PLAN-serve-index-and-cleanup.md`. Differences from this draft: the UI starts from one query object and four actions (annotation TSV, loci TSV/BED, GFF slice link, FASTA); `transcript` resolution moves server-side; saved gene sets, "Use in API" and the compare/heatmap extras are optional later items.

Scope: `annotrieve/front/app/genes/page.tsx` and `annotrieve-xrefs/api/`. v1 is gene symbols only. The end goal is the FASTA (gffread) service.

This is an opinion based on reading the page, the `/hits`, `/names`, `/annotations/...` and `/sequences` routers, `helpers/concurrency.py`, `helpers/pagination.py`, the start of `services/sequences.py`, the pin store and `lib/genes/sequences-download.ts`. It has not been tested. The rest of `services/sequences.py` (including the Ensembl accession mapping) and the Compare view were not read.

## 1. UI flow

**The tabs are organized by view, but the end goal is a task.** Today a user gets sequences by searching, going to Annotations, opening the drawer, pinning a gene, switching to Compare, waiting for transcript resolution, and then clicking Download. That is about six hops, capped at 40 pins. FASTA is the point of the feature and it is buried.

Ideas:
- **Search -> Scope -> Select -> Export as a visible flow.** Keep the three tabs as ways to look at the results. Add a persistent selection tray and a primary "Get sequences" button that is available in every tab.
- **Two selection modes.**
  - *Explicit pins*, as now.
  - *Everything in the current filtered view*, such as "all 143 annotations matching BRCA1, one per species, RefSeq preferred". The user should not have to pin 143 rows. This needs the API change below.
- **Export presets** for the common cases: CDS / protein / transcript, longest isoform vs all isoforms, one annotation per species. `keep.ts` may already do the "one per species" part. If it does, reuse it for export.
- **A preflight dialog before the job starts.** For example: "187 genes, 143 species, about 40 MB, about 6 min, 12 will be skipped (reason)". Then show progress, then the download plus a skip report.
- **Honest coverage messaging.** Absence of a symbol is not absence of the gene. Names differ across taxa (LOC ids, lineage-specific naming, paralogs). The zero-hit review is a good start. A line like "found in 62% of annotations; unmatched ones may use a different name" would stop users over-interpreting.
- **Flag ambiguous matches.** Matching is symbol ∪ alias, so an alias can hit a paralog. Show whether each hit matched by symbol or by alias, and flag annotations with several loci for the same query. `locus-siblings.ts` suggests part of this is handled. Check that it is visible at export time.
- **Other entry points.**
  - A "find gene in this annotation" link from the annotation and organism pages.
  - A "Use in API" box with copyable curl and Python for the current query.
  - Short saved gene sets (`POST /sets` returns an id), because long lists do not fit in URLs and reproducible exports need a stable id.

## 2. API shape and logic

**Biggest problem: `/sequences` is one long synchronous request.**
- The deadline is 800 s, there is no progress reporting, and any failure loses everything (fail-first). Reverse proxies and browsers will time out on long exports.
- Proposed shape:
  - `POST /sequence-jobs` returns `202 {job_id}`.
  - `GET /sequence-jobs/{id}` returns the phase (slice, NCBI fetch, gffread), `done/total`, and the skip list.
  - `GET /sequence-jobs/{id}/result` returns the FASTA.
- Keep the sync path for small requests (for example 10 items or fewer).
- Cache results by a hash of (sorted items, action, buffer).
- Add `on_error: "fail" | "skip"`, with skipped items listed in a manifest. One unresolvable seqid in 200 items should not kill the export.

**Selection is item-only, which forces client-side N+1 calls.**
- The client needs a `local_id` per gene via `/annotations/{id}/genes`, then a resolved `transcript_id`, before it can call `/sequences`.
- Let `/sequences` accept a query instead: `{names, annotation_ids | filter, per_species: "one", transcript: "longest" | "all" | <id>}`. The server expands it using the same logic as `/hits`.
- Also allow `transcript: "longest"` on explicit items. The store already treats `null` as "use longest", but the API requires an id.

**Add `POST /sequences/estimate` (dry run).** It returns the gene count, the number of bp to fetch from NCBI, the cache hits, the unresolvable items, and an estimated time. It powers the preflight dialog and lets you refuse oversized jobs early.

**The NCBI rate limit is probably not enforced globally.**
- `concurrency.py` uses `threading.Semaphore`, which is per gunicorn worker. With W workers you get W x `SEQUENCES_MAX_INFLIGHT` x `NCBI_WORKERS` concurrent efetches against one NCBI key.
- You need a shared token bucket (Redis, or a file lock/sqlite) or a single dedicated sequence worker. Do this before launch, with or without the higher NCBI limit.
- The slice cache is only 32 MB in memory. A persistent on-disk cache keyed by `accession:start-end` would cut NCBI traffic a lot.
- Long term, a local mirror of the genomes users request most would remove the dependency.

**Output options:**
- `grouping: "single" | "by_gene_zip"`. For alignments and phylogenetics, people want one file per gene with one sequence per species.
- A defline template or presets that include species, symbol, annotation id and assembly.
- A TSV manifest alongside the FASTA, with one row per sequence.

**Smaller items:**
- `/hits` could optionally return per-term gene copy counts, so the UI can show "BRCA1 x2" before export.
- The per-annotation shard routes (xrefs, namespaces) are not used by the symbol-only v1. Leave them but do not advertise them.
- The 200-item cap is fine for explicit-item mode. With query-based selection, enforce limits on the server-side expansion instead (total bp or sequence count).
- Add a per-IP quota on the sequence endpoints, since they hit a shared external resource.

## 3. Suggested order
1. Shared NCBI rate limiter plus persistent slice cache. This is correctness and a prerequisite for everything else.
2. `estimate` plus the job model with `on_error: skip` and a manifest.
3. Query-based selection with server-side `longest` and one-per-species.
4. A UI "Get sequences" CTA with preflight and progress, replacing the pin -> Compare -> Download path (keep Compare for visualization).
5. Saved gene sets and the "Use in API" box.

## 4. Open questions
- Do Ensembl-sourced annotations always map to NCBI accessions (`resolve_seqid_accession`)? If not, those exports will fail, and the UI should say so before the job starts.
- What is the largest realistic export: a few genes across hundreds of species, or hundreds of genes in a few species? The first favours per-species dedup and by-gene zips. The second favours a local genome mirror.
- Is `/sequences` meant to be public and scripted, or UI-only for now? That decides how much quota and auth work is needed.
