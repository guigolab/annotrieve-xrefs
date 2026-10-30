# annotrieve-xrefs

Search [Annotrieve](https://genome.crg.es/annotrieve) annotations by **gene name** (symbol / synonym), then drill into genes and xrefs inside one annotation.

**Base URL:** `https://genome.crg.es/annotrieve-xrefs/api/v1`

This service is a sidecar to Annotrieve. Annotrieve owns the main annotation catalogue; this API owns **name → annotation** presence search and per-annotation gene/xref drill-in. It is not a replacement for NCBI Gene or Ensembl.

Typical flow (happy path):

1. **Discover names (optional)** — `GET /names?prefix=…` (prefix suggestions with categories)
2. **Search by exact gene names** — `GET|POST /hits?names=…` (symbol∪alias; no wildcards)
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


More than **500** distinct names → `422 prefix_too_broad`. No pagination (bounded by that expand cap). Results ordered by `n_annotations` desc, then `name`. Categories (first-match-wins): `pseudogene`, `numbered`, `homeolog`, `predicted`, `exact`, `other`.

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

**Gene-name presence search** (primary API). Exact names only against **symbol ∪ alias** (casefolded). A name containing `*` → `422 wildcard_not_supported` (use `GET /names` for prefix discovery).


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

### `POST /sequences`

Bulk transcript FASTA export (POST only; never cached). For each item the API resolves the gene from the shard, finds the requested transcript in the on-disk GFF (one tabix query per gene), fetches only that transcript’s genomic window from NCBI E-utilities, and runs **gffread** to emit CDS, spliced transcript, or protein sequence. Fail-first: the first unresolved seqid / missing transcript / NCBI / gffread error aborts with **no partial body**.

```http
POST /sequences
Content-Type: application/json

{
  "action": "cds",
  "buffer_bp": 0,
  "include_exon_coords": false,
  "items": [
    {
      "annotation_id": "ann-human",
      "local_id": 10,
      "transcript_id": "XM_015139305.2"
    }
  ]
}
```

| Field | Notes |
| --- | --- |
| `action` | Required. `cds` (gffread `-x`), `transcript` (spliced exons+UTRs, `-w`), or `protein` (`-y`). |
| `items` | Required. 1..**200** `{annotation_id, local_id, transcript_id}`. `transcript_id` is always required (no longest-transcript default). Duplicates are dropped, first-seen order kept. |
| `buffer_bp` | Optional. 0..**10000** (default 0). Extra flank on each side; passed to gffread as `--w-add`. Only valid with `action=transcript` (`422 buffer_not_supported` otherwise). |
| `include_exon_coords` | Optional bool (default false). Passes gffread `-W`. |

**200** `text/x-fasta` with `Content-Disposition: attachment` and `X-Sequence-Count`. Deflines look like:

```text
>XM_015139305.2|TP53|Homo_sapiens|ann-huma NC_000017.11:7661779-7687550(-) buffer=0
```

Items are grouped by gene: several isoforms of one gene cost one tabix query, one NCBI fetch (union of requested spans ± buffer), and one gffread run. Overlapping windows on the same accession within 10 kb are coalesced.

**Errors** use the shared `{"detail": {"code", "message", …}}` shape (extras include `index`, `annotation_id`, `local_id`, `transcript_id` when relevant):

| Status | Codes |
| --- | --- |
| 422 | `validation_error`, `too_many_items`, `items_required`, `buffer_not_supported`, `seqid_not_resolved`, `no_cds`, `transcript_incomplete`, `transcript_outside_gene`, `transcript_strand_unknown`, `window_too_large` |
| 404 | `shard_not_found`, `gene_not_found`, `annotation_not_found`, `annotation_file_not_found`, `transcript_not_found` |
| 503 | `sequences_busy` (+ `Retry-After: 1`), `annotations_unavailable`, `ncbi_rate_limited` (+ `Retry-After`), `gffread_unavailable` |
| 502/504 | `ncbi_sequence_not_found`, `ncbi_error`, `ncbi_timeout`, `ncbi_sequence_truncated`, `sequences_timeout` |
| 500 | `gffread_failed` |

**Ops / env** (see `.env.example`; never put the NCBI key in compose):

- `ANNOTATIONS_DIR` — mount of Annotrieve’s annotation files (prod: `/data/annotations`).
- `ANNOTATIONS_GFF_PATH_PREFIX` — optional harvest `--files-root` to rebase `annotation.gff_path`.
- `NCBI_API_KEY` — optional; with key ~9 req/s total, without ~2.5 req/s (split across `NCBI_WORKERS`).
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

