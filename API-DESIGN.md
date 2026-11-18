# Annotrieve xrefs: API design (v1, symbol and alias only)

Status: **implemented and live-tested** (2026-11): `/names`, `/hits`, `/loci` and `/sequences` (with `/sequence-estimates`, `/sequence-plans`, `/sequence-jobs`) all read the presence and loci indexes (`PRESENCE_DIR`; `loci.sqlite` built by `serve_build loci`, `loci_ref.sqlite` by `loci-ref`). Shared across them: `alias=off|suggest|include`, taxon scope, prefix terms (`seleno*`, with `exclude`, `categories`, `split`), strict parameters (unknown parameter -> `422 unknown_parameter`). `/hits` has `sources`, `min_matched`, `shape=compact|records`, `format=json|tsv` and `facets`; `/names` has `sources` and `n_loci`; `/loci` has json, tsv, bed and ndjson with a completeness trailer. Still open: ETag/cache headers, the sync mechanism (D6), namespaces (D7), UI changes. This is the reference the UI and the build pipeline follow. Companion documents:
- `PLAN-build-tables.md`: the script and logic that build the new tables (sidecars with coordinates, a manifest, `presence.sqlite`, `loci.sqlite`). Built and run on the real corpus.
- `PLAN-serve-index-and-cleanup.md`: the served files, memory budget, and what to drop.
- `REVISED_PLAN.md`: earlier design (organisms, genome grain). Parts are superseded here.

Sequencing (historical): presence first, then `/loci` and `/sequences`; both are done. The sync mechanism for new, replaced and removed annotations (D6) is still to be decided.

## 1. What the service is for

It closes the gap between **genes** and **GFF files** in Annotrieve. A user has gene names and wants:

1. **Annotations:** which of the ~18k annotations contain them, as a table or TSV with download links.
2. **Coordinates:** where the matching genes are, so they can fetch GFF slices from Annotrieve.
3. **Sequences:** FASTA (CDS, transcript or protein) of those genes via gffread.

Out of scope for v1: non-symbol namespaces (GO, Pfam, InterPro, ...), orthologs, free-text search. They stay in the shards but are not served.

## 2. Principles

1. **One query object drives everything.** The same terms, scope and filters are used by `/names`, `/hits`, `/loci` and `/sequences`. The UI keeps one state, and every export is reproducible from it.
2. **Exports are driven by the query, not by id lists.** The server expands the query. Explicit id lists exist only for small selections.
3. **Formats are a parameter:** `format=json` (default), `tsv`, and `bed` where it makes sense. TSV and BED start with `#` provenance lines (query, expanded terms, `data_version`, timestamp).
4. **Coordinates:** JSON and TSV use 1-based inclusive (GFF convention). BED uses 0-based half-open. Every locus also carries a ready `gff_url`, so users never convert by hand. Annotrieve's slice endpoint (`GET /annotations/{md5}/gff?region=&start=&end=`) documents `start` as 0-based and `end` as inclusive, so `gff_url` must be built and tested against it. The service does not re-implement slicing.
5. **Every response carries `data_version`,** used as the ETag.
6. **Refuse before doing the work.** Every list reports `total`. An oversized export returns `413` with an estimate, computed from counts and never from a scan.
7. **GET is canonical and cacheable** (lowercase, sorted, deduplicated names). POST mirrors GET for long lists.
8. **Stable references.** Annotations are identified by `annotation_id` (md5 of the file, stable). A gene is referenced as `locus_ref = annotation_id + ":" + gene key`. See the open decision on `gene_id` in section 9.

## 3. The query object and prefix expansion

```json
{
  "terms": [
    {"q": "tp53"},
    {"q": "seleno*", "categories": ["exact"], "exclude": ["selenoh"], "split": false}
  ],
  "include_taxids": [40674],
  "exclude_taxids": [],
  "sources": ["RefSeq", "Ensembl"],
  "min_matched": 2
}
```

GET shorthand: `names=tp53,seleno*&include_taxids=40674&sources=RefSeq,Ensembl&min_matched=2`.

### Terms
- A **term** is either an **exact name** (`tp53`) or a **prefix** (`seleno*`, a single trailing `*`, at least 3 characters).
- A prefix is expanded **by the server** against the sorted name dictionary (a range scan, never `LIKE`).
- A term matches an annotation if **any** of its expanded names (symbol or alias) is present. Coverage (`n_matched`, `min_matched`) counts **terms**, not names. This keeps "I asked for 2 things" meaningful when one of them is a prefix.
- `split: true` turns each expanded name into its own term, so coverage is counted per name.
- Per-term annotation counts are exact **unions** of the names' presence, never sums.

### Matching semantics: symbol, alias, case
- **A match means the gene's symbol.** `n_matched`, `min_matched`, coverage and the TSV counts use symbol hits only.
- **Alias hits are a suggestion, not a match.** Aliases are only about 2.8% of the index rows, but they carry the cross-species false positives (for example the alias `oor` belongs to `Oprl1` in some annotations). For each term the response reports how many **more** annotations match only through an alias (`n_alias_extra`), and `alias=suggest` (default) returns those annotations flagged `alias_only`, with the gene's own primary name. Parameter `alias=off|suggest|include`: `off` ignores aliases, `suggest` (default) keeps them out of coverage but lists them, `include` restores the old symbol-or-alias union.
- **Case.** Names are indexed casefolded, so `tp53`, `TP53` and `Tp53` find the same annotations by default. Original-case spellings come from the gene's primary name (about 95% of symbol rows). `case_sensitive=true` is applied to `/loci` rows (exact match on the primary name) and is not available on `/hits` presence in v1; whether to index exact spellings is decided from the spelling check in the plan.

### Prefix options
- `categories`: restrict the expansion to name categories (`exact`, `numbered`, `homeolog`, `pseudogene`, `predicted`, `other`), as the prefix sheet does today.
- `exclude`: names to drop from the expansion (the user unticking rows in the sheet).
- Expansion is scoped: the counts shown for each name are for the current taxon scope.

### Limits
- Prefix: at least 3 characters. At most 500 names per prefix term. At most 2,000 expanded names per request (POST).
- Over the limit: `422 prefix_too_broad` with `{n_names, n_annotations, max}`. These numbers come from the name dictionary alone, so a refused request is cheap (today it costs seconds).

### Reproducibility
Responses and export headers echo the expansion:

```json
"terms": [
  {"index": 0, "input": "tp53", "kind": "exact", "names": ["tp53"], "n_annotations": 8440, "n_alias_extra": 12},
  {"index": 1, "input": "seleno*", "kind": "prefix", "names": ["selenop", "selenow", "..."],
   "n_names": 37, "n_annotations": 3120, "truncated": false}
]
```

### UI flow stays the same
`/names` is still the discovery call behind the prefix sheet (scoped counts and categories per name). The user picks names, and the committed set can be sent either as exact names or as the prefix term plus its `exclude` list. The second is shorter and keeps the intent ("everything starting with seleno") visible.

## 4. Endpoints

### 4.1 `GET /names?prefix=&include_taxids=&exclude_taxids=&sources=` (discovery)
Distinct names under a prefix that are a **symbol** somewhere in scope (alias-only names, mostly locus-style IDs, are listed only with `alias=include`), with `category` and `n_annotations` in scope. **Implemented.** `n_loci` per name (gene rows in the same scope; symbol rows only unless `alias=include`; absent when `loci.sqlite` is not published), `sources=Ensembl,GenBank,RefSeq` to count only those annotations (echoed as `source_filter`), and `422 prefix_too_broad` above 500 names.

### 4.2 `GET|POST /hits` (which annotations) - **implemented**

Parameters: the query object (`names` and, on POST, `terms`), `include_taxids`, `exclude_taxids`, `sources`, `min_matched`, `alias=off|suggest|include` (default `suggest`), `shape=compact|records` (default `compact`, the UI's positional `fields`/`rows`), `format=json|tsv` (`tsv` implies records and streams everything). No `limit`/`offset`: the whole list is returned (sorted by terms matched, then id).

Record (JSON):

```json
{
  "annotation_id": "3f1c...",
  "organism_name": "Homo sapiens", "taxid": 9606,
  "assembly_accession": "GCF_000001405.40", "assembly_name": "GRCh38.p14",
  "source": "RefSeq", "release_date": "2024-04-01",
  "n_matched": 2, "n_terms": 2, "matched": [0, 1], "alias_only": [], "missing": [],
  "download_url": "..."
}
```

`shape=compact` keeps today's positional `fields` and `rows` for the UI, now with the added columns.

TSV columns: `annotation_id, organism_name, taxid, assembly_accession, assembly_name, source, release_date, n_matched, n_terms, matched_terms, alias_only_terms, missing_terms, download_url`.

Response also carries: `terms` (expansion echo), `total`, `facets` (counts by source), `data_version`.

### 4.3 `GET|POST /loci` (where the genes are) - **implemented** (deviations: records list `terms`, `matched_names` and `matched_via` for the whole gene instead of a single `term`/`matched_name`; a `gene_id` column holds the source gene id; `limit` on tsv/bed/ndjson stops the stream after N rows, and without it unscoped queries above `LOCI_MAX_ROWS` get 413; `locus_ref` = `annotation_id:local_id`; `gff_url` uses pysam's 0-based half-open slice, so GFF S..E becomes `start=S-1&end=E`)

Parameters: the query object, plus `annotation_ids` (optional restriction), `flank_bp` (applied to `gff_url`), `format=json|tsv|bed`, `dry_run=true`, `alias=off|suggest|include` (default `suggest`), `case_sensitive` (default false), `limit` and `cursor` for JSON.

Record (JSON):

```json
{
  "annotation_id": "3f1c...", "locus_ref": "3f1c...:123456",
  "gene_id": null, "symbol": "TP53",
  "term": 0, "matched_name": "tp53", "matched_via": "symbol",
  "seqid": "NC_000017.11", "start": 7661779, "end": 7687550, "strand": "-",
  "feature_type": "gene", "biotype": "protein_coding",
  "gff_url": "https://genome.crg.es/annotrieve/api/v0/annotations/3f1c.../gff?region=NC_000017.11&start=7661778&end=7687550"
}
```

TSV columns: `annotation_id, organism_name, assembly_accession, source, symbol, term, matched_name, matched_via, seqid, start, end, strand, feature_type, biotype, locus_ref, gff_url`.
BED columns: `seqid, start-1, end, name = symbol|annotation_id, score = 0, strand`. Note: `seqid` is the annotation's own name (accession or chromosome label). The Annotrieve slice endpoint accepts aliases.

`dry_run=true` returns `{n_annotations, n_loci_upper_bound}` and no rows.

Pipelines want the whole result set, so exports are streamed, not paginated:
- `format=tsv|bed|ndjson` stream everything in deterministic order (annotation ordinal, then seqid, then start). Rows are produced in bounded batches, so memory does not grow with the result.
- A final trailer line (`#end rows=N complete=true`, or an `error` object in NDJSON) lets a client detect a truncated stream.
- A hard cap (default 2M rows, a placeholder until measured) returns `413` with the estimate before streaming starts; the client narrows the query or raises the limit with `limit=`.
- JSON (UI use) is paginated with `limit` and `cursor`.

One locus per gene: a gene matching several terms appears once with the matched terms listed.

### 4.4 `POST /sequences` (FASTA) - **implemented** (new index only, no old shard or corpus table; needs `loci_ref.sqlite` for `loci` lists, built by `serve_build loci-ref`)

Flat JSON body. Exactly one source of genes:

```json
{ "query": { "names": ["tp53"], "include_taxids": [40674], "alias": "suggest" },
  "action": "cds", "transcript": "longest", "on_error": "skip" }
```
```json
{ "loci": ["2fcd45eb...:1043", {"locus_ref": "a1b2...:77", "transcript_id": "ENST00000269305"}],
  "action": "protein" }
```

| field | meaning |
|---|---|
| `query` | same fields as `/loci`; reproducible from the UI state |
| `loci` | list of `locus_ref` strings, or objects with a pinned `transcript_id` (skips the rule). Replaces the legacy `items` (kept one release) |
| `action` | `cds` \| `transcript` \| `protein` (required) |
| `transcript` | preset name or rule object (below). Default `longest` |
| `buffer_bp`, `include_exon_coords` | as today; `buffer_bp` only with `action=transcript` |
| `on_error` | `skip` (default for `query`) or `fail`; skipped genes go to the manifest |
| `max_genes` | cap for a `query`; over the cap -> refused with the estimate, never silently cut |

**Transcript choice** (`helpers/transcript_select.py`, implemented and tested). There is deliberately **no canonical** option: the flag exists in one data source only (measured 2026-11-18: Ensembl 96/96 files; RefSeq only in human and mouse; GenBank and community none).

Presets: `longest` (default), `longest_spliced`, `longest_cds`, `shortest_cds`, `most_exons`, `fewest_exons`, `all`.
`longest` depends on the action: `cds`/`protein` -> longest CDS; `transcript` -> coding first, then longest spliced length. `longest_spliced` is the literal longest, coding or not.

Rule object: `{"rank_by": ["coding:desc", "cds_len:desc"], "require": ["coding", "tag:basic", "biotype:protein_coding"], "keep": 1}`.
- Rank keys: `coding`, `cds_len`, `spliced_len`, `span`, `n_exons`, `utr_len`; each `:asc` or `:desc` (bare = desc).
- `require`: `coding`, `tag:<value>`, `biotype:<value>`. `keep`: 1..50 or `"all"`.
- Eligibility comes first: for `cds`/`protein` only transcripts with CDS are considered; for `transcript` any with structure.
- Ties always go to the lowest transcript id (not configurable), so a request is reproducible.
- Invalid input: `unknown_transcript_preset`, `unknown_rank_key`, `unknown_require`, `invalid_transcript_rule` (each lists the allowed values).

**Implementation notes** (deviations and choices, 2026-11-18)
- Batches of `SEQUENCES_BATCH_GENES` (25) genes bound memory to one batch of NCBI windows; progress is reported per batch. Results live in `SEQUENCES_RESULTS_DIR/<result_id>/` (`result.fa`, `manifest.tsv`, `status.json`), so any worker serves them and a job is polled the same way from any worker; they expire after `SEQUENCES_RESULT_TTL_S`.
- `job_id` equals `result_id`. Up to `SEQUENCES_SYNC_MAX_GENES` (50) genes run in the request; more become a job (`202`, `Location`). `max_genes` (default and ceiling `SEQUENCES_MAX_GENES`, 2000) is checked by counting the genes first: over it is `413 too_many_genes`.
- `on_error` defaults to `skip` for `query` and `fail` for `loci`/`items`. Skip reasons also include `annotation_not_found`, `locus_not_found`, `gene_id_missing`, `transcript_strand_unknown`, `gffread_failed`, `ncbi_sequence_truncated`.
- A transcript that reaches beyond its gene bounds makes the resolver fetch the wider GFF region again, so exons are never silently missing (the old pipeline failed such genes with `transcript_outside_gene`).
- `locus_ref` lookup needs `loci_ref.sqlite` (`ord, local_id -> name_id`, one row per gene, built with `serve_build loci-ref` or automatically by `serve_build loci`); without it, `loci` lists return `503 locus_ref_index_unavailable` and `query` mode still works.
- FASTA header: `>{symbol}|{organism}|{assembly}|{annotation_id}|{transcript_id} {seqid}:{start}-{end}({strand})`, with spaces and `|` inside a field replaced by `_`.
- Done since: the NCBI limiter shared by all workers and a per-IP quota (below). Not done: `transcript` ids typed as a plain string (`transcript=ENST...`; use a pinned `loci` entry instead).

**Companion resources** (stateless, computed from the same body as `POST /sequences`; nothing is stored; renamed from `/sequences/estimate` and `/sequences/plan`; typed in `api/openapi.json`)
- `POST /sequence-estimates`: index only, no GFF reads: `{n_genes, n_annotations, est_bp_upper_bound, over_limit}`. Transcript counts are unknown at this stage.
- `POST /sequence-plans`: same body, reads the GFF slices. One row per eligible transcript with its metrics, `selected`, and the rule that applied. A client edits the table and posts it back as `loci`.

**Storage and abuse limits.** Results and scratch files live under `SEQUENCES_RESULTS_DIR` (shared storage, TTL 24 h), at most `SEQUENCES_MAX_RESULTS` (100) at a time. Each result records a hash of the client IP; a client may have 2 exports running and start 20 per hour (429 `too_many_active_exports` / `export_rate_limited`, with `Retry-After`). The check is soft (two simultaneous requests can both pass) and is computed from the retained results, so every worker agrees. The NCBI rate-limit file stays on local `/tmp` because `flock` is unreliable on NFS.

**NCBI budget.** All API processes share one request budget through a lock file (`NCBI_RATE_STATE_FILE`); a 429 from NCBI delays every process. Per-client fairness is not enforced beyond the proxy's per-IP limit.

**Response.** Every result has a `result_id` = hash of the normalised request + `data_version`.
- Small requests (default up to 50 genes): `200 text/x-fasta` with headers `X-Sequence-Count`, `X-Skipped-Count`, `X-Result-Id`, `X-Data-Version`.
- Larger requests: `202 {job_id, result_id}`, then `GET /sequence-jobs/{id}` (phase, `done/total`).
- `GET /sequences/results/{result_id}` (FASTA) and `/manifest.tsv` while the result is cached. FASTA has no safe place for comments, so skip details live only in the manifest.
- Manifest columns: `locus_ref, symbol, status, transcript_id, rule, n_candidates, reason, seq_length`.
- Skip reasons: `no_transcript`, `no_cds`, `filtered_out`, `transcript_not_found`, `gff_unavailable`, `seqid_not_resolved`, `window_too_large`, `ncbi_sequence_not_found` (the contig is not an NCBI sequence), `ncbi_timeout`, `ncbi_error`.
- FASTA header: `>{symbol}|{organism}|{assembly_accession}|{annotation_id}|{transcript_id}`.
- Transcript models come from the GFF slice when a request needs them (no global transcript table), with a read-through cache keyed by `annotation_id` and gene. GFF path comes from the presence index (`annotation.taxid`, `annotation.gff_path`).
- Limits: see "Storage and abuse limits" and "NCBI budget" above (shared NCBI limiter, per-IP quota, results cached by `result_id`).

### 4.5 Kept as they are
`GET /annotations/{id}/genes`, `/genes/{local_id}`, `/xrefs`, `/namespaces`: single-annotation drill-in. They read the shards, so they stay for gene cards and are not used for bulk work.

## 5. Errors, caching, limits

- Error body `{"detail": {"code", "message", ...}}` as today. New codes: `prefix_too_broad` (now with counts), `export_too_large` (413 with estimate), `job_not_found`.
- Today `data_version` is in the body and in an `X-Data-Version` header on streamed responses. **Not done:** `ETag` from `data_version` and a canonical nginx cache key.
- Concurrency gates and `503 Retry-After` stay.
- Rate limit at the edge as now.

## 6. What each endpoint reads

| Endpoint | Reads |
|---|---|
| `/names` | `name_dict` range scan; postings only for scoped counts of at most 500 names |
| `/hits` | `name_dict`, `posting`, in-memory `annotation` and lineage |
| `/loci` | `name_dict` (for expansion), then `locus` range scans per name; `annotation` for labels |
| `/sequences` | `locus` (to expand the query), then the existing tabix, gffread and NCBI pipeline |

Table definitions and sizes are in `PLAN-serve-index-and-cleanup.md`.

## 7. What the UI follows

The page already has three views (`components/genes/`): **Distribution** (heatmap of genes across taxa), **Annotations** (list and table) and **Compare** (a genome-browser-style view that draws the transcript structures of pinned genes side by side, with a toolbar that already downloads CDS, transcript or protein FASTA through `postSequences`). Only minor changes are needed; the views stay.

- One query bar produces the query object and the same object feeds all three views.
- **Annotations:** `/hits?format=tsv` is the annotation-link TSV; `/loci?format=tsv` is the one with the locus and `gff_url`. The UI has not been switched to them yet.
- **Compare / FASTA:** keep the pin tray and transcript picker; add the estimate before a large request and a job/progress state, and a "longest vs all transcripts" choice.
- **Aliases and prefixes** are visible and reversible in the UI (alias suggestions as a separate strip; prefix expansion preview with unticking), never silent.
- Freshness (indexed N of M, pending) appears once the sync mechanism is decided.

## 8. Migration from the current API

- `GET /hits` keeps working. `shape=compact` is today's response plus the new columns; `names=` without `*` behaves as now.
- A `*` in `names` is a prefix term (done; it was `wildcard_not_supported`). Without the presence index (legacy mode) `*` is still rejected.
- `/sequences` still accepts `items` with `annotation_id, local_id, transcript_id` for one release, then moves to `locus_ref`.

## 9. Open decisions

1. **Gene identifier in loci (decided: `locus_ref = annotation_id:local_id`, rows also carry `gene_id`).** The fast source for coordinates, `gene_hit`, has `local_id` but not the source gene id (`gene.source_gene_id` lives in the shards). Options: (a) `locus_ref` = annotation id + `local_id` with the caveat that a rebuilt shard can change it; (b) read `source_gene_id` from the shards during the build (slow: the shards total 418 GiB); (c) fill `gene_id` lazily from the shard when a user opens one gene or runs a FASTA job. Recommended for v1: (a) plus (c), and (b) only for new annotations as they arrive.
2. **Annotation download URL.** Decided: `download_url` is `{ANNOTRIEVE_FILES_BASE}/{taxid}/{assembly}/{file}` (the last three parts of the stored path). It goes stale when an assembly's taxid folder changes until the next index sync (D6). `gff_url` offsets are tested against the slice endpoint.
3. **Loci file size.** Measured on the real index: 80,995,350 genes (see `PLAN-build-tables.md`, loci build).
4. **Deleted, replaced and new annotations.** Search must exclude removed ones at once and include new ones after the weekly import. The mechanism (a reconcile job against the live Annotrieve catalogue, a liveness mask, a `GET /status`) is **postponed** until the first table works; see `PLAN-build-tables.md` D6.
5. **Sync size cut-offs** for `/sequences` (50 genes) and jobs (2,000 genes). Placeholders until measured.
6. **"Longest transcript" rule**: decided 2026-11-18, see 4.4 (action-aware, ranking keys, no canonical option). Open: where transcript models are cached.
7. **Namespaces and the annotation filter** (Pfam, GO, ...) are deferred; see `PLAN-build-tables.md` D7. Only Ensembl, GenBank and RefSeq annotations are indexed; `Community` is skipped.
