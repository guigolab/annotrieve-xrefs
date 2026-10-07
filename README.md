# annotrieve-xrefs

Search Annotrieve annotation files by gene symbol, gene id, or other accession, and see the matching genes.

Public API: `https://genome.crg.es/annotrieve-xrefs/api/v1`

This service is a sidecar to [Annotrieve](https://genome.crg.es/annotrieve). It serves a published SQLite index (`gene_corpus.sqlite` + per-annotation `genes.sqlite` shards) with keyset pagination (`next` / `previous`). Harvest and merge live in `gene_corpus/`.

## Facets

| Facet | Endpoint |
|---|---|
| Discover | `GET /namespaces`, `GET /namespaces/{namespace}/accessions` |
| Quick search | `GET /genes?id=symbol:tp53&id=alias:tp53&taxid=40674` |
| Details | `GET /annotations/{annotation_id}/genes/{local_id}` |

## Local development

```bash
# Optional: point at a published corpus directory
export GENE_CORPUS_HOST_DIR=/path/to/gene_corpus

docker compose -f docker-compose-dev.yml up --build
```

- API: http://localhost:5003/health
- Via nginx: http://localhost:95/annotrieve-xrefs/api/v1/health

## Harvest / merge

```bash
python -m gene_corpus.harvest --work-dir DIR --files-root FILES
python -m gene_corpus.merge --work-dir DIR \
  --taxonomy-tsv /path/to/flattened-tree.tsv
```

`--taxonomy-tsv` is Annotrieve’s flattened taxonomy export (`taxid`,
`parent_taxid` columns; e.g. `/annotrieve/files/taxonomy/flattened-tree.tsv`).
Merge streams `gene_hit`, upserts `xref_meta` (`n_annotations`, `n_loci`), and
fills `annotation_lineage`. Peak SQLite heap is soft-capped at 512 MiB.
Interrupted merges can continue with `--resume` / `--resume-after-key`.

## Tests

```bash
pip install -r requirements-dev.txt
PYTHONPATH=api:. pytest gene_corpus/tests api/tests -q
```
