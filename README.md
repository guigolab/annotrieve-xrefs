# annotrieve-xrefs

Search [Annotrieve](https://genome.crg.es/annotrieve) annotations by cross-reference — gene symbol, Gene Ontology, Pfam, UniProt, and related feature ids.

**Base URL:** `https://genome.crg.es/annotrieve-xrefs/api/v1`

This service is a sidecar to Annotrieve. Annotrieve owns the main annotation catalogue; this API owns **feature → annotation** search and per-annotation gene/xref drill-in. It is not a replacement for NCBI Gene or Ensembl.

Typical flow:

1. **Discover** — `GET /matches?q=` (which namespaces match) then `GET /namespaces/{ns}/terms?q=` (pick CURIEs)
2. **Search annotations** — `GET|POST /annotations?curies=…` (presence; compact heatmap rows)
3. **Drill in** — `GET /annotations/{id}/genes…` and `…/xrefs`

---

## Pagination

Paged list endpoints use **`limit`**, **`offset`**, and exact **`total`**:

```json
{ "total": 814, "limit": 100, "offset": 0, "results": [ … ] }
```

- Default `limit` is **100**; maximum is **200**.
- Maximum `offset` is **10_000**.

Bulk annotation catalogue / presence search use a compact `fields` + `rows` payload (no offset paging).

Failed responses: `{"detail": {"code": "...", "message": "..."}}`. Validation failures are **422** with `code=validation_error`.

Production rate limits: **3 req/s per IP**, ~**20 req/s** globally. Saturated annotation search concurrency returns **503** with `Retry-After: 1`.

---

## CURIEs

Search ids use Compact URIs: `prefix:accession` (first colon splits prefix from accession).

| Prefix in the request | Stored namespace |
|---|---|
| `swiss-prot` | `uniprot/swiss-prot` |
| any searchable name (e.g. `symbol`, `go`, `pfam`) | same name |

Normalization:

- **GO** — short forms pad to seven digits (`GO:8150` → `GO:0008150`).
- **symbol / alias** — accessions are casefolded (`TP53` → `tp53`).

Only **searchable** namespaces are accepted as CURIE prefixes. Unknown prefixes appear in the `errors` array on `/annotations?curies=…` and on `/annotations/…/genes?curies=…` (valid CURIEs still return under `match=any`).

---

## Searchable namespaces

API catalogue and CURIE search are limited to:

`symbol`, `alias`, `geneid`, `go`, `pfam`, `interpro`, `rfam`, `mirbase`, `cdd`, `eggnog`, `uniprot/swiss-prot`, `ncbiortholog`

Other harvested namespaces (Ensembl, model-organism DBs, TrEMBL, …) may still appear on **per-gene xref chips** inside a shard, but are not global search facets.

Namespaces with a label collection (GO, Pfam, InterPro, Rfam, miRBase, Swiss-Prot) report `"labels": true` on `/namespaces` and `/matches`.

---

## API reference

Paths below are relative to the base URL.

### `GET /namespaces`

Directory of searchable namespaces present in the corpus. Fixed shape (no `q`, no `match_count`).

| Parameter | Notes |
|---|---|
| `limit`, `offset` | Pagination |

```json
{
  "total": 4,
  "limit": 100,
  "offset": 0,
  "results": [
    {"namespace": "alias", "name": "Gene synonym", "accession_count": 1, "labels": false},
    {"namespace": "go", "name": "Gene Ontology", "accession_count": 1, "labels": true},
    {"namespace": "symbol", "name": "Gene symbol", "accession_count": 1, "labels": false}
  ]
}
```

### `GET /namespaces/{namespace}`

One namespace summary (same fields as a directory row). Slash namespaces work in the path (e.g. `/namespaces/uniprot/swiss-prot`).

### `GET /matches`

**Namespace facet** for a query string. `q` is required. Every result includes **`match_count`**.

Labeled namespaces match on label text **or** accession (Pfam free-text accepts versioned or unversioned ids).

| Parameter | Notes |
|---|---|
| `q` | CURIE, numeric id, symbol token, or prose (min 3 chars except CURIE / digits) |
| `limit`, `offset` | Pagination over matching namespaces |

```http
GET /matches?q=kinase
```

```json
{
  "total": 2,
  "limit": 100,
  "offset": 0,
  "results": [
    {"namespace": "go", "name": "Gene Ontology", "accession_count": 23999, "labels": true, "match_count": 512},
    {"namespace": "interpro", "name": "InterPro", "accession_count": 60651, "labels": true, "match_count": 814}
  ]
}
```

### `GET /namespaces/{namespace}/terms`

Paged terms in one searchable namespace.

| Parameter | Notes |
|---|---|
| `q` | Optional for labeled namespaces (browse; FTS or accession); required for symbol/alias/geneid/… |
| `limit`, `offset` | Pagination; `total` is exact |

```json
{
  "namespace": "go",
  "total": 1,
  "limit": 100,
  "offset": 0,
  "results": [
    {
      "curie": "go:GO:0008150",
      "namespace": "go",
      "accession": "GO:0008150",
      "label": "biological_process",
      "n_annotations": 2,
      "n_genes": 2
    }
  ]
}
```

### `GET /namespaces/{namespace}/terms/{accession}`

One term card (`curie`, `label`, counts). GO / symbol forms are normalized like CURIEs.

### `GET` / `POST /annotations`

Without `curies`: **placement catalogue** for the heatmap (annotations with at least one Tier A xref in `xref_file`).

With `curies`: **presence search** — annotations that carry the CURIEs (`match=any|all`, optional `taxid`). Same compact payload.

```http
GET /annotations
GET /annotations?curies=symbol:tp53&taxid=9606&match=any
POST /annotations
{"curies": ["symbol:tp53", "GO:0008150"], "match": "all", "taxid": 9606}
```

```json
{
  "sources": {"1": "Ensembl", "2": "GenBank", "3": "RefSeq"},
  "fields": ["annotation_id", "annotation_key", "taxid", "source"],
  "rows": [["ann-human", 1, 9606, 3]]
}
```

Catalogue-only response also includes `organism_name` and `assembly_accession` in `fields` / `rows`, plus `n_without_tier_a`.

Join presence rows to the catalogue on `annotation_key`. Built by `python -m gene_corpus.index`.

### `GET /annotations/{annotation_id}/namespaces`

Namespaces present on one annotation’s shard (Tier A and chip namespaces), with distinct accession counts.

### `GET /annotations/{annotation_id}/genes`

Browse / filter genes inside one annotation.

| Parameter | Notes |
|---|---|
| `q` | Casefold **prefix** on `primary_name` (max 64 chars) |
| `curies` | Optional CSV CURIEs via `gene_xref` (cap 20) |
| `match` | `any` (default) or `all` (with `curies`) |
| `limit`, `offset` | Pagination; response includes exact `total` |

### `GET /annotations/{annotation_id}/genes/{local_id}`

One gene card (coordinates, name, prose).

### `GET /annotations/{annotation_id}/genes/{local_id}/xrefs`

Paginated xref chips for one gene (`limit` / `offset` / `total`).

---

## Building the corpus

1. **Harvest** — per-annotation gene / xref shards  
2. **Merge** — global `gene_corpus.sqlite` (catalogue, meta, optional `gene_hit`)  
3. **Index** — `xref_file` presence for `/annotations`  
4. **Labels** (optional) — `xref_label` / FTS for prose discovery  

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

- API: http://localhost:5003
- Via nginx: http://localhost:95/annotrieve-xrefs/api/v1

### Tests

```bash
pip install -r requirements-dev.txt
PYTHONPATH=api:. pytest gene_corpus/tests api/tests -q
```
