# annotrieve-xrefs

Search [Annotrieve](https://genome.crg.es/annotrieve) annotations by **gene name** (symbol / synonym), then drill into genes and xrefs inside one annotation.

**Base URL:** `https://genome.crg.es/annotrieve-xrefs/api/v1`

This service is a sidecar to Annotrieve. Annotrieve owns the main annotation catalogue; this API owns **name → annotation** presence search and per-annotation gene/xref drill-in. It is not a replacement for NCBI Gene or Ensembl.

Typical flow (happy path):

1. **Search by gene names** — `GET|POST /hits?names=…` (symbol∪alias; optional `*` prefix)
2. **Drill in** — `GET /annotations/{id}/genes?names=…` and `…/xrefs`
3. **Or browse xrefs** — `…/namespaces` → `…/namespaces/{ns}/accessions?q=…` → `…/genes?curies=ns:acc`

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

### `GET` / `POST /hits`

**Gene-name presence search** (primary API). Names are matched against **symbol ∪ alias** (casefolded). Optional trailing `*` expands a prefix (min 3 characters before `*`; at most 500 accessions or `422 prefix_too_broad`).


| Parameterel      | Notes                                                                                                                   |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `names`          | Required. CSV (GET) or JSON list (POST). Cap GET **50** / POST **200**.                                                 |
| `include_taxids` | Optional CSV/list; keep annotations whose lineage includes any listed taxid. Cap **50** taxids (`422 too_many_taxids`). |
| `exclude_taxids` | Optional; drop annotations whose lineage includes any listed taxid. Cap **50**.                                         |
| `match`          | `any` (default) or `all` — union vs intersection over **input names**                                                   |


```http
GET /hits?names=SELENOP,SEPP1,gpx4*&include_taxids=7742&match=any
POST /hits
{"names": ["SELENOP", "SEPP1", "gpx4*"], "match": "all", "include_taxids": [9606]}
```

```json
{
  "terms": [
    {"input": "SELENOP", "name": "selenop", "matched_namespaces": ["symbol"], "n_hits": 1},
    {"input": "gpx4*", "name": "gpx4*", "matched_namespaces": ["symbol"],
     "expanded": ["gpx4", "gpx4a", "gpx4b"], "n_hits": 2}
  ],
  "keys": [
    {"term": 0, "name": "selenop", "ns": "symbol"},
    {"term": 1, "name": "gpx4", "ns": "symbol"}
  ],
  "sources": {"1": {"label": "Ensembl"}, "2": {"label": "GenBank"}, "3": {"label": "RefSeq"}},
  "include_taxids": [7742],
  "exclude_taxids": [],
  "total": 2,
  "fields": ["annotation_id", "taxid", "source", "n_matched", "matched"],
  "rows": [
    ["ann-human", 9606, 3, 2, [0, 1]]
  ]
}
```

`matched` lists indices into `keys`. `n_matched` is the number of distinct input **terms** those keys cover. `matched_namespaces` / `n_hits` on each term reflect postings **after** tax filters. The response is the full hit set (no pagination); clients derive taxon sets from row `taxid`s.

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

Browse / filter genes inside one annotation. Use **`names`** (symbol∪alias; same token rules as `/hits`, optional `*` prefix) when drilling in from `/hits`, or **`curies`** (`ns:accession` exact keys for any shard namespace) after picking from `/namespaces/{ns}/accessions`. `names` and `curies` are mutually exclusive (`422 names_and_curies`). Taxid filters stay on `/hits` only. `annotation_id` is path-only (not repeated in the JSON).


| Parameter         | Notes                                                                                          |
| ----------------- | ---------------------------------------------------------------------------------------------- |
| `names`           | Optional CSV gene names (cap **50**). Omit to browse all genes. Mutually exclusive with `curies`. |
| `curies`          | Optional CSV CURIEs (cap **20**). Exact accessions only. Mutually exclusive with `names`.       |
| `match`           | `any` (default) or `all` — union vs intersection over **input names or CURIEs**                 |
| `limit`, `offset` | Pagination; response includes exact `total`                                                    |


When `names` is set, each gene includes `matches`: `{input, name, ns}`. When `curies` is set: `{input, ns, accession}`. Omitted when browsing without a filter.

```http
GET /annotations/ann-human/genes?names=SELENOP,gpx4*&match=any&limit=50
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
export GENE_CORPUS_HOST_DIR=/path/to/gene_corpus
docker compose -f docker-compose-dev.yml up --build
```

- API: [http://localhost:5003](http://localhost:5003)
- Via nginx: [http://localhost:95/annotrieve-xrefs/api/v1](http://localhost:95/annotrieve-xrefs/api/v1)



### Tests

```bash
pip install -r requirements-dev.txt
PYTHONPATH=api:. pytest gene_corpus/tests api/tests -q
```

