# gene_corpus

The offline half of annotrieve-xrefs. It reads annotation files (GFF3), keeps the gene names and
cross-references, and builds the small SQLite files that the API (`../api`) serves. Nothing here
runs inside the API process; it runs as two command-line steps on the host, with flat memory
(under about 2 GB) and scratch files on the NFS work area.

```
Annotrieve report ──harvest──▶ per_gff/<id>/genes.sqlite      (one shard per annotation)
                                    │
                    serve_build catalogue / sidecars
                                    ▼
              per_gff/<id>/sidecar.tsv.gz  +  serve/manifest.sqlite
                                    │
                    serve_build build / loci / loci-ref
                                    ▼
        serve/v<data_version>/{presence,loci,loci_ref}.sqlite   ◀── serve/current (symlink)
```

## The two commands

| Command | What it does |
|---|---|
| `python -m gene_corpus.harvest` | Reads the Annotrieve annotation report, opens each GFF3 (local file or URL), and writes one shard `per_gff/<annotation_id>/genes.sqlite` (genes, their cross-references, per-namespace counts, metadata). Options: `--work-dir`, `--gff-mode local\|url`, `--files-root`, `--annotation-id` / `--annotation-ids-file` (only those), `--workers`, `--memory-flush-units`, `--force`, `--api-base`, `--report`. |
| `python -m gene_corpus.serve_build <sub-command>` | Turns shards into the served files. Sub-commands, in order: `catalogue` (load annotation metadata into `serve/manifest.sqlite`), `sidecars` (one sorted `sidecar.tsv.gz` per shard), `build` (`presence.sqlite`, published atomically as `serve/current`), `loci` (`loci.sqlite` next to it), `loci-ref` (`loci_ref.sqlite`, lookup by `locus_ref`). Helpers: `missing` (catalogue annotations with no shard yet) and `status` (manifest counts). |

Both take `--work-dir` (the corpus root that holds `per_gff/` and `serve/`).

## Top-level files

| File | Role |
|---|---|
| `__init__.py` | Package docstring and the few names re-exported (`canonicalize_db`, `include_in_local`, `include_in_v1`). |
| `harvest_cli.py` | Implements `python -m gene_corpus.harvest`: reads the report, picks jobs, runs workers, calls `stream.harvest` per annotation. |
| `serve_build_cli.py` | Thin dispatcher for `python -m gene_corpus.serve_build`; it only forwards to `serve_backfill_cli`. |
| `serve_backfill_cli.py` | The real implementation of every `serve_build` sub-command (argument parsing, workers, publishing). The name is historical. |
| `gating.py` | Per-file parse gates: from the report's metadata (attribute keys, feature types, roots) decide what a GFF is worth parsing. |
| `resolve.py` | Finds the GFF (local path or URL) for an annotation and picks the right profile (RefSeq, GenBank or Ensembl). |

## Folders

| Folder | Role |
|---|---|
| `harvest/` | Entry point only: `__main__.py` calls `harvest_cli`, so `python -m gene_corpus.harvest` works. |
| `serve_build/` | Entry point only: `__main__.py` calls `serve_build_cli`. |
| `catalog/` | The annotation report from Annotrieve: `report.py` downloads and parses it (`POST /annotations/report`, filtered to Ensembl, GenBank and RefSeq); `models.py` holds `AnnotationJob` and `ParseGate`. |
| `namespaces/` | Which cross-reference databases count and under what name. `aliases.py` canonicalises Dbxref/source labels and flags noise (specimen, taxonomy); `tiers.py` defines the tier A/B namespace sets (also used by the API for `kind`/`slug`); `labels.py` maps a declared label to a namespace. |
| `profiles/` | One class per source (`refseq.py`, `genbank.py`, `ensembl.py`) that says which GFF column-9 attributes become symbols, aliases, ids and cross-references, and which are ignored. |
| `stream/` | Reading and parsing a GFF3 in one pass. `reader.py` streams plain or gzip lines and skips unwanted landmarks early; `attrs.py` parses column 9; `prose.py` picks the gene description; `emit_util.py` shared emit helpers (GO ids, Ensembl descriptions); `symbol_gate.py` decides whether a name is a real symbol or alias; `harvest.py` writes one shard, flushing open loci in bounded windows. |
| `sql/` | Everything that writes SQLite or sidecar files (see below). |
| `tests/` | The unit tests of the folders above (and of `tools/free_disk.py`). Run with `python3 -m pytest gene_corpus`. |

### `sql/`

| File | Role |
|---|---|
| `per_gff.py` | Schema of a shard: tables `gene`, `gene_xref`, `tier_a_counts`, `meta`. Written once; replaced as a whole on rebuild. |
| `sidecar.py` | Writes the per-annotation `sidecar.tsv.gz`: sorted gene names (symbol or alias) and the named genes with their coordinates. Fingerprints tell when a sidecar is stale. |
| `manifest.py` | `manifest.sqlite`: which annotations are indexed, their stable ordinals (never reused), status and counts. `work_shard` runs in workers; `record_result` runs in the single process that owns the manifest. |
| `catalogue.py` | Loads Annotrieve's annotation metadata (organism, taxid, database, release date, lineage) into the manifest, streaming line by line. |
| `postings.py` | Codec for posting lists (which annotations contain a name): varint deltas or bitset, whichever is smaller. |
| `serve_presence.py` | Shared pieces of the served-index builds: the `presence.sqlite` schema, `BuildError`, and the atomic `publish_version` switch. |
| `presence_from_sidecars.py` | Builds `presence.sqlite` (name to annotations, plus lineage) by an external merge of the sidecars, in flat memory. |
| `loci_from_sidecars.py` | Builds `loci.sqlite` (one row per name and gene, with coordinates) and `loci_ref.sqlite` (lookup by `locus_ref`), and checks them against the presence file name by name. |

## Files the build produces

```
<work-dir>/per_gff/<annotation_id>/genes.sqlite       harvest
<work-dir>/per_gff/<annotation_id>/sidecar.tsv.gz     serve_build sidecars
<work-dir>/serve/manifest.sqlite                      serve_build catalogue / sidecars
<work-dir>/serve/v<data_version>/presence.sqlite      serve_build build
<work-dir>/serve/v<data_version>/loci.sqlite          serve_build loci
<work-dir>/serve/v<data_version>/loci_ref.sqlite      serve_build loci-ref
<work-dir>/serve/current -> v<data_version>           switched atomically by build
```

The API opens `serve/current` and reopens it when the link changes. To remove old versions, partial
builds and other leftovers, use `python3 tools/free_disk.py` (dry run by default).
