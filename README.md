# annotrieve-xrefs

Search [Annotrieve](https://genome.crg.es/annotrieve) annotations by **gene name** (symbol / synonym), then drill into genes and xrefs inside one annotation.

**Base URL:** `https://genome.crg.es/annotrieve-xrefs/api/v1`

This service is a sidecar to Annotrieve. Annotrieve owns the main annotation catalogue; this API owns **name → annotation** presence search and per-annotation gene/xref drill-in. It is not a replacement for NCBI Gene or Ensembl.

Typical flow (happy path):

1. **Discover names (optional)** — `GET /names?prefix=…` (prefix suggestions with categories)
2. **Search by gene names** — `GET|POST /hits?names=…` (symbol∪alias; a trailing `*` makes a prefix term)
3. **Drill in** — `GET /annotations/{id}/genes?names=…` and `…/xrefs`
4. **Or browse xrefs** — `…/namespaces` → `…/namespaces/{ns}/accessions?q=…` → `…/genes?curies=ns:acc`

---



## Pagination

Paged list endpoints use `limit`, `offset`, and exact `total`:

```json
{ "total": 814, "limit": 100, "offset": 0, "results": [ … ] }
```

- Default `limit` is **100**; maximum is **200**.
- Maximum `offset` is **10_000**.

`/hits` returns the full hit set (compact `fields` + `rows`; no offset paging).

Failed responses: `{"detail": {"code": "...", "message": "..."}}`. Validation failures are **422** with `code=validation_error`.

Production rate limits: **3 req/s per IP**, ~**20 req/s** globally. Saturated `/hits` concurrency returns **503** with `Retry-After: 1`.

---



## Namespaces

`/hits` and `/annotations/{id}/genes?names=…` match **symbol ∪ alias** (casefolded; `TP53` → `tp53`).

Inside one annotation, browse any harvested namespace (Tier A and chip namespaces) via `/namespaces` → `/namespaces/{ns}/accessions`, then resolve with `/genes?curies=ns:accession`. Per-gene chips remain at `…/genes/{local_id}/xrefs`.

---



## API reference

Paths bow are relative to the base URL.

### `GET /names`

**Prefix discovery** for the gene-name sheet. Returns distinct symbol∪alias accessions under a prefix, with annotation counts in the current taxon scope and a server-assigned `category`.


| Parameter        | Notes                                                                                                                   |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `prefix`         | Required. Casefolded stem, no `*`. At least **3** characters (`422 prefix_too_short`). `*` → `422 invalid_prefix`.     |
| `include_taxids` | Optional CSV; keep annotations whose lineage includes any listed taxid. Cap **50** (`422 too_many_taxids`).             |
| `exclude_taxids` | Optional; drop annotations whose lineage includes any listed taxid. Cap **50**.                                         |


`sources` (Ensembl, GenBank, RefSeq; comma-separated) counts only those annotations. With the loci index published each result also has `n_loci` (gene rows in the same scope; symbol rows only unless `alias=include`), so a client can warn before a large `/loci` export; without it the field is absent. More than **500** distinct names → `422 prefix_too_broad`. No pagination (bounded by that expand cap). Results ordered by `n_annotations` desc, then `name`. Categories (first-match-wins): `pseudogene`, `numbered`, `homeolog`, `predicted`, `exact`, `other`.

```http
GET /names?prefix=seleno&include_taxids=7742
```

```json
{
  "prefix": "seleno",
  "include_taxids": [7742],
  "exclude_taxids": [],
  "total": 37,
  "results": [
    {"name": "selenop", "category": "exact", "n_annotations": 8123}
  ]
}
```

### `GET` / `POST /hits`

**Gene-name presence search** (primary API). Names are matched against **symbol ∪ alias** (casefolded). A term ending in `*` (`seleno*`, same rule on `/hits`, `/loci` and the `query` of `/sequences`) is a **prefix term**: the server expands it against the name dictionary (at least 3 characters; at most 500 names per term and 2,000 per request, else `422 prefix_too_broad`; alias-only names are included only with `alias=include`) and echoes the expansion in `terms` (`kind`, `names`, `n_names`). A term matches if any of its names does, and `min_matched` counts terms. POST accepts `terms`: strings or objects `{q, exclude, categories, split}` (options for prefix terms only; `split` turns each expanded name into its own term; `categories` as in `GET /names`). `*` elsewhere → `422 invalid_prefix`.


| Parameter        | Notes                                                                                                                   |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `names`          | Required. CSV (GET) or JSON list (POST). Cap GET **200** / POST **200**. Casefolded duplicates count once.             |
| `include_taxids` | Optional CSV/list; keep annotations whose lineage includes any listed taxid. Cap **50** taxids (`422 too_many_taxids`). |
| `exclude_taxids` | Optional; drop annotations whose lineage includes any listed taxid. Cap **50**.                                         |


```http
GET /hits?names=SELENOP,SEPP1,TP53&include_taxids=7742
POST /hits
{"names": ["SELENOP", "SEPP1", "TP53"], "include_taxids": [9606]}
```

```json
{
  "terms": [
    {"input": "SELENOP", "name": "selenop", "namespaces": ["symbol"], "n_hits": 1},
    {"input": "TP53", "name": "tp53", "namespaces": ["symbol"], "n_hits": 1}
  ],
  "sources": {"1": {"label": "Ensembl"}, "2": {"label": "GenBank"}, "3": {"label": "RefSeq"}},
  "include_taxids": [7742],
  "exclude_taxids": [],
  "total": 2,
  "fields": ["annotation_id", "taxid", "source", "matched", "release_date", "assembly_name"],
  "rows": [
    ["ann-human", 9606, 3, [0, 1], "2024-01-02T00:00:00", "GRCh38.p14"]
  ]
}
```

`matched` lists **term indices** into `terms`. `namespaces` / `n_hits` on each term reflect postings **after** tax filters. The response is the full hit set (no pagination); clients derive taxon sets from row `taxid`s. Prefer GET with a canonical query string so nginx can cache identical sets.

### `GET /annotations/{annotation_id}/namespaces`

Namespaces present on one annotation’s shard (Tier A and chip namespaces), with distinct accession counts. `annotation_id` is path-only (not repeated in the JSON).

### `GET /annotations/{annotation_id}/namespaces/{namespace}/accessions`

Browse / prefix-search accessions under one namespace on this shard. Path `{namespace}` accepts a stored namespace, a public short prefix (`swiss-prot` → `uniprot/swiss-prot`), or a slash form (`uniprot/swiss-prot` or `uniprot%2Fswiss-prot`). For `symbol` / `alias`, `q` is casefolded to match stored accessions. **404** `namespace_not_found` when the namespace is absent on the shard.


| Parameter         | Notes                                                                                          |
| ----------------- | ---------------------------------------------------------------------------------------------- |
| `q`               | Optional accession **prefix**. Omit to browse all. Non-empty must be ≥ **3** characters (`422`). |
| `limit`, `offset` | Pagination; response includes exact `total`                                                    |


Each result: `{accession, display, gene_count}` (`display` nullable; representative via `MAX(display)` when several loci share an accession).

```http
GET /annotations/ann-human/namespaces/refseq_transcript/accessions?q=XM_&limit=50
```

```json
{
  "total": 1,
  "limit": 50,
  "offset": 0,
  "results": [
    {
      "accession": "XM_015139305.2",
      "display": "XM_015139305.2",
      "gene_count": 1
    }
  ]
}
```

### `GET /annotations/{annotation_id}/genes`

Browse / filter genes inside one annotation. Use **`names`** (exact symbol∪alias; wildcards rejected) when drilling in from `/hits`, or **`curies`** (`ns:accession` exact keys for any shard namespace) after picking from `/namespaces/{ns}/accessions`. `names` and `curies` are mutually exclusive (`422 names_and_curies`). Taxid filters stay on `/hits` only. `annotation_id` is path-only (not repeated in the JSON).


| Parameter         | Notes                                                                                            |
| ----------------- | ------------------------------------------------------------------------------------------------ |
| `names`           | Optional CSV exact gene names (cap **200**). Omit to browse all genes. Mutually exclusive with `curies`. |
| `curies`          | Optional CSV CURIEs (cap **20**). Exact accessions only. Mutually exclusive with `names`.         |
| `match`           | `any` (default) or `all` — union vs intersection over **input names or CURIEs**                   |
| `limit`, `offset` | Pagination; response includes exact `total`                                                      |


When `names` is set, each gene includes `matches`: `{input, name, ns}`. When `curies` is set: `{input, ns, accession}`. Omitted when browsing without a filter.

```http
GET /annotations/ann-human/genes?names=SELENOP,TP53&match=any&limit=50
GET /annotations/ann-human/genes?curies=refseq_transcript:XM_015139305.2
```

```json
{
  "total": 1,
  "limit": 50,
  "offset": 0,
  "results": [
    {
      "local_id": 10,
      "primary_name": "SELENOP",
      "matches": [
        {"input": "SELENOP", "name": "selenop", "ns": "symbol"}
      ]
    }
  ]
}
```



### `GET /annotations/{annotation_id}/genes/{local_id}`

One gene card (coordinates, name, prose). No `annotation_id` field (path-scoped).

### `GET /annotations/{annotation_id}/genes/{local_id}/xrefs`

Paginated xref chips for one gene (`limit` / `offset` / `total`). `annotation_id` and `local_id` are path-only (not repeated in the JSON).

### `/downloads` (FASTA, queued)

Ask for up to **1,000** CDS, transcript or protein sequences. The work runs in **one queue for the whole service** (one download at a time, because of the NCBI request budget), so a request never waits on the connection. The body is the one described in `API-DESIGN.md` 4.4 (`query` or `loci`, `action`, `transcript`, `buffer_bp`, `on_error`, ...).

| Step | Request | Answer |
| --- | --- | --- |
| Check the size first | `POST /downloads?dry_run=true` | `200`: gene count, base-pair bound, limits. Nothing is queued and no quota is used. |
| Ask | `POST /downloads` | `202` with `Location` and the status. An equal request returns the download that exists (`200` when finished) and costs no quota. |
| Follow | `GET /downloads/{id}` | `state` (`queued`, `running`, `done`, `failed`), `queue_position`, `progress`, `links`; `Retry-After` while waiting |
| Fetch | `GET /downloads/{id}/fasta` | gzip-compressed FASTA (`application/gzip`); `409 download_not_ready` before `done` |
| Skipped genes | `GET /downloads/{id}/report` | tab-separated outcome per gene |
| Cancel or delete | `DELETE /downloads/{id}` | `204`; only the client that asked |
| Your quota | `GET /downloads/quota` | `limit`, `used`, `remaining`, `resets_at`, `active` |

Limits per client IP: 10 new downloads per rolling 24 hours (`429 daily_quota_exceeded`, with `Retry-After`) and 2 waiting or running at once (`429 too_many_active_downloads`). Files are removed 10 minutes after the first complete download, 2 hours after `done` if never fetched, 1 hour after a failure. The store holds at most 100 downloads and the queue 50 (`503 downloads_queue_full` / `downloads_capacity`). FASTA deflines: `>{symbol}|{organism}|{assembly}|{annotation_id}|{transcript_id} {seqid}:{start}-{end}({strand})`.

**Ops / env** (see `.env.example`; never put the NCBI key in compose):

- `ANNOTATIONS_DIR` — mount of Annotrieve’s annotation files (prod: `/data/annotations`).
- `PRESENCE_DIR` — directory with the `current` symlink to `v<version>/presence.sqlite` (built by `serve_build build`). When set, `/names` and `/hits` read it and support `alias=off|suggest|include` (default `suggest`: symbol hits are matches, alias-only annotations are listed as `alias_only`, `n_alias_extra` per term/name). Unset = legacy `xref_file` behaviour. `/loci` also needs `loci.sqlite` in the same version directory; related settings: `ANNOTRIEVE_API_BASE` (for `gff_url`), `LOCI_MAX_ROWS`, `LOCI_CACHE_KIB`, `LOCI_MAX_INFLIGHT`. `/downloads` also needs `ANNOTATIONS_DIR` (GFF files), `gffread` and, for `loci` lists, `loci_ref.sqlite` in the version directory (`serve_build loci-ref`); settings: `SEQUENCES_MAX_GENES` (1000) and `SEQUENCES_MAX_SEQUENCES` (1000): the most one download holds (over it: `413 too_many_genes` or `413 too_many_sequences`, before queueing), `SEQUENCES_DEADLINE_S` (1800: longest a download may run), `SEQUENCES_PLAN_MAX_GENES` (200), `SEQUENCES_BATCH_GENES` (25), `SEQUENCES_RESULTS_DIR` (downloads, quota ledger and scratch, shared by all workers; in compose the host folder `/data/annotrieve-xrefs/sequences` mounted at `/data/results`, which must be writable by uid 1000), `DOWNLOADS_PER_DAY` (10 new downloads per client IP per rolling 24 h; equal requests are free; 0 disables), `SEQUENCES_IP_MAX_ACTIVE` (2 waiting or running per client IP), `DOWNLOADS_MAX_QUEUED` (50), `SEQUENCES_MAX_RESULTS` (100 stored; nothing is removed to make room, new downloads get 503), `DOWNLOAD_KEEP_AFTER_FETCH_S` (600), `DOWNLOAD_KEEP_UNFETCHED_S` (7200), `DOWNLOAD_KEEP_FAILED_S` (3600), `DOWNLOADS_LOCK_FILE` (the one-runner lock, local disk: `/tmp/annotrieve-downloads.lock`), `DOWNLOADS_POLL_S` (2), `TRUSTED_PROXY_HOPS` (1: the client IP is the last `X-Forwarded-For` entry, so the nginx template must keep passing that header); `NCBI_RATE_STATE_FILE` (all API processes share one NCBI budget, 2.5 req/s or 9 with `NCBI_API_KEY`, through a lock file; default `<SEQUENCES_RESULTS_DIR>/.ncbi-rate`; a 429 from NCBI slows every process; if the file cannot be opened each process limits itself to rate / `NCBI_PROCESSES`, set that to the gunicorn worker count); `PUBLIC_BASE_URL` (public scheme, host and path prefix, e.g. `https://genome.crg.es/annotrieve-xrefs/api/v1`) so `Location` and the links are correct behind the proxy.
- `ANNOTATIONS_GFF_PATH_PREFIX` — optional harvest `--files-root` to rebase `annotation.gff_path`.
- `NCBI_API_KEY` — optional; with key ~9 req/s, without ~2.5 req/s (one shared token bucket per API process).
- `NCBI_WORKERS` — fetch thread-pool size only (does not divide the rate).
- `NCBI_EMAIL`, `NCBI_TOOL`, gunicorn `--timeout 900`, nginx `proxy_read_timeout 900s` on this location.

---



## Building the corpus

1. **Harvest** — per-annotation gene / xref shards
2. **Merge** — global `gene_corpus.sqlite` (catalogue, meta, optional `gene_hit`)
3. **Index** — `xref_file` presence for `/hits`
4. **Labels** (optional) — `xref_label` / FTS (corpus tooling; not exposed as HTTP discovery)

```bash
PYTHONPATH=. python -m gene_corpus.index --work-dir /path/to/gene_corpus
PYTHONPATH=. python -m gene_corpus.labels --work-dir /path/to/gene_corpus
```

---



## Appendix: local development

```bash
cp .env.example .env   # optional NCBI_API_KEY=
export GENE_CORPUS_HOST_DIR=/path/to/gene_corpus
export ANNOTATIONS_HOST_DIR=../annotrieve/files
docker compose -f docker-compose-dev.yml up --build
```

- API: [http://localhost:5003](http://localhost:5003)
- Via nginx: [http://localhost:95/annotrieve-xrefs/api/v1](http://localhost:95/annotrieve-xrefs/api/v1)



### Tests

```bash
pip install -r requirements-dev.txt
PYTHONPATH=api:. pytest gene_corpus/tests api/tests -q
```

## API schema

The OpenAPI schema is committed as `api/openapi.json` (interactive docs at `/docs` when the service runs). After changing a route or a model, regenerate it with `python3 tools/export_openapi.py`; a test fails when the file is stale. Unknown top-level fields in a `/sequences` body are rejected with `validation_error` instead of being ignored.
