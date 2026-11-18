# Annotrieve xrefs: API design (v1, symbol and alias only)

Status: **proposal; the presence build (script and tests) and the presence-backed `/names` and `/hits` with `alias=off|suggest|include` exist** (**`/loci` exists too**: exact names, `alias`, taxon scope, `annotation_ids`, `flank_bp`, `case_sensitive`, `format=json|tsv|bed|ndjson` with a completeness trailer, `dry_run`, JSON cursor, 413 export cap; it needs `loci.sqlite`, built by `serve_build loci`) (compact `fields`/`rows` kept as the only shape for now, with `alias_only`, `organism_name`, `assembly_accession` appended; enabled by `PRESENCE_DIR`; prefix terms in `/hits`, `shape`, `format=tsv`, `min_matched` and `facets` not done yet). This is the reference the UI and the build pipeline follow. Companion documents:
- `PLAN-build-tables.md`: the script and logic that build the new tables (sidecars with coordinates, a manifest, `presence.sqlite`). **Current focus.**
- `PLAN-serve-index-and-cleanup.md`: the served files, memory budget, and what to drop.
- `REVISED_PLAN.md`: earlier design (organisms, genome grain). Parts are superseded here.

Sequencing: the first table (`presence.sqlite`, behind `/names` and `/hits`) comes first. `/loci` (its input, coordinates, is already in the sidecar) and the sync mechanism for new, replaced and removed annotations are decided afterwards. Sections below that depend on them are marked **later**.

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
Distinct names under a prefix that are a **symbol** somewhere in scope (alias-only names, mostly locus-style IDs, are listed only with `alias=include`), with `category` and `n_annotations` in scope. Added: `n_loci` per name (so the UI can warn before a large loci export), and the refusal path described above.

### 4.2 `GET|POST /hits` (which annotations)

Parameters: the query object, plus `alias=off|suggest|include` (default `suggest`), `shape=records|compact` (default `records`), `format=json|tsv`, `limit`/`offset` (JSON records only; TSV streams everything).

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

### 4.4 `POST /sequences` (FASTA) - **later** (needs `/loci`); the current endpoint keeps working meanwhile

```json
{
  "query": { "...": "same query object as above" },
  "items": null,
  "action": "cds",
  "transcript": "longest",
  "buffer_bp": 0,
  "include_exon_coords": false,
  "on_error": "skip"
}
```

- Exactly one of `query` (all matching loci) or `items` (explicit `locus_ref` list, small selections).
- `transcript`: `longest` (default, resolved on the server when the job reads the GFF lines), `all`, or a transcript id (with `items`).
- `action`: `cds | transcript | protein`; `buffer_bp` only with `transcript`.
- `on_error`: `fail` (today's behaviour) or `skip`. Skipped items go to a manifest.
- `POST /sequences/estimate` returns `{n_genes, est_bp_upper_bound, est_seconds, over_limit}`, computed from the loci table (gene length is an upper bound for the bp to fetch).
- Small requests (default up to 50 genes) return the FASTA directly with `X-Sequence-Count` and `X-Skipped-Count`. Larger requests return `202 {job_id}`, then `GET /sequence-jobs/{id}` (phase and `done/total`) and `GET /sequence-jobs/{id}/result` (FASTA) and `/manifest.tsv`.
- FASTA header: `>{symbol}|{organism}|{assembly_accession}|{annotation_id}|{transcript_id}`.
- Transcripts are **not precomputed** in a global table (the corpus is large and keeps growing). Transcript models come from the GFF slice when a request needs them (as the Compare view does today), with a read-through cache keyed by `annotation_id` and gene. "Longest" is resolved from those models; the exact length rule (genomic span, summed exons or CDS length) is still to decide.
- Limits: a global NCBI rate limiter shared across all workers, a per-IP quota, and results cached by a hash of the normalized request plus `data_version`.

### 4.5 Kept as they are
`GET /annotations/{id}/genes`, `/genes/{local_id}`, `/xrefs`, `/namespaces`: single-annotation drill-in. They read the shards, so they stay for gene cards and are not used for bulk work.

## 5. Errors, caching, limits

- Error body `{"detail": {"code", "message", ...}}` as today. New codes: `prefix_too_broad` (now with counts), `export_too_large` (413 with estimate), `job_not_found`.
- nginx cache key is the canonical query string. `data_version` appears in `ETag`.
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
- **Annotations:** the export becomes a TSV with the locus and a way to fetch the chunk (annotation id plus `gff_url`), once `/loci` exists. Until then it stays the annotation-link TSV from `/hits`.
- **Compare / FASTA:** keep the pin tray and transcript picker; add the estimate before a large request and a job/progress state, and a "longest vs all transcripts" choice.
- **Aliases and prefixes** are visible and reversible in the UI (alias suggestions as a separate strip; prefix expansion preview with unticking), never silent.
- Freshness (indexed N of M, pending) appears once the sync mechanism is decided.

## 8. Migration from the current API

- `GET /hits` keeps working. `shape=compact` is today's response plus the new columns; `names=` without `*` behaves as now.
- A `*` in `names` stops being an error and becomes a prefix term.
- `/sequences` keeps accepting `items` with `annotation_id, local_id, transcript_id` for one release, then moves to `locus_ref`.

## 9. Open decisions

1. **Gene identifier in loci.** The fast source for coordinates, `gene_hit`, has `local_id` but not the source gene id (`gene.source_gene_id` lives in the shards). Options: (a) `locus_ref` = annotation id + `local_id` with the caveat that a rebuilt shard can change it; (b) read `source_gene_id` from the shards during the build (slow: the shards total 418 GiB); (c) fill `gene_id` lazily from the shard when a user opens one gene or runs a FASTA job. Recommended for v1: (a) plus (c), and (b) only for new annotations as they arrive.
2. **Annotation download URL.** Confirm the file-link pattern Annotrieve uses for a whole annotation, and the exact offset in the slice endpoint's `start`.
3. **Loci file size.** A guess of 2.5-3.5 GB. Measure by building it, and time a cold range read for a common name first.
4. **Deleted, replaced and new annotations.** Search must exclude removed ones at once and include new ones after the weekly import. The mechanism (a reconcile job against the live Annotrieve catalogue, a liveness mask, a `GET /status`) is **postponed** until the first table works; see `PLAN-build-tables.md` D6.
5. **Sync size cut-offs** for `/sequences` (50 genes) and jobs (2,000 genes). Placeholders until measured.
6. **"Longest transcript" rule** (genomic span, summed exons, or CDS length) and where transcript models are cached.
7. **Namespaces and the annotation filter** (Pfam, GO, ...) are deferred; see `PLAN-build-tables.md` D7. Only Ensembl, GenBank and RefSeq annotations are indexed; `Community` is skipped.
