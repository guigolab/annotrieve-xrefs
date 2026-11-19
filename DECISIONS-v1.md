# V1 decisions

Status: draft for V1 (2026-11-19). This is the scope contract for the xrefs API and the genes page. Each decision has a **trigger**: the evidence that reopens it. Until a trigger fires, the answer to "can it do X?" is "not in V1".

## The job

> Give a list of gene names (symbols, aliases, prefixes) and a taxonomic scope. Get back where those genes are, in which annotations, and their sequences, in bulk, with provenance.

The gap this fills: NCBI Datasets and Ensembl BioMart answer "gene X in species Y" well, one source at a time. Nothing offers a uniform bulk lookup by gene name across thousands of eukaryotic annotations from every source (Ensembl, RefSeq, GenBank, community). The value is **breadth and uniformity**, not depth per gene.

## Rules for where to stop

1. **Critical path.** A feature is V1 only if the job needs it (name, then presence, then loci, then FASTA, with provenance) or if it makes results honest.
2. **Global index or shard.** Index globally only what a query needs across annotations. Everything else stays in the per-annotation shards and is served as detail on demand. Keeping harvested data in shards is cheap; a new global index is the expensive decision.
3. **No competing with upstream.** GO enrichment, orthology inference, domain search and sequence search belong to other tools. Link to them.
4. **Evidence before indexes.** The API has no users yet. Ship, log real queries, then build where the logs and users point.

## In V1

- `/names`, `/hits`, `/loci` (json, tsv, ndjson) and `/sequences` with `/sequence-estimates`, `/sequence-plans`, `/sequence-jobs`.
- Symbols and aliases (gated by `gene_corpus/stream/symbol_gate.py`), prefix terms, taxon include and exclude, `alias=off|suggest|include`.
- Per-annotation gene lists and gene detail, including xrefs, read from the shards.
- Provenance on every output: `data_version`, annotation, source, release.
- A visible **coverage statement**: how many annotations are searchable by name, per source (for example "9,545 of 17,199"). Without it, "not found" is misleading.

## Deferred, with triggers

| Feature | Decision | Revisit when |
|---|---|---|
| Search by GO, Pfam or InterPro | Defer. Data stays in shards, no global index. | Users report nameless annotations as their blocker. Strongest v2 candidate: domains are how a family is found in annotations without names. |
| Lookup by Ensembl or GeneID | Defer. It would be a point-lookup table (id to annotation and gene), not a presence index. | ID input shows up in the logs or in user requests. |
| NCBI ortholog identifiers and synonym groups | Defer. Partial coverage, and easy to misread as orthology. | A user needs name-independent cross-species grouping. |
| Transcript table | Start with a small per-gene summary (transcript count, longest CDS, exon count, flags), not the full table. | The UI needs bulk gene structure for thousands of rows. |
| Species-list scope | Likely V1.1. | First user feedback. |
| Relaxing the name gate | No. Decision of 2026-11-18 stays. | Never, unless re-measured with `tools/measure_gate_relax.py`. |

## Measured facts behind these decisions

From `tools/measure_corpus_size.py` on the 2026-11-18 build (`data_version` 9e8dc2f914f42ba5):

- 17,199 annotations (GenBank 9,494, Ensembl 5,110, RefSeq 2,595) over 9,002 taxa.
- 9,545 annotations (55%) have indexed names. The rest have no gate-passing symbol or alias, mostly locus-style names (`LOC123`, `TE_00001670`, ab initio IDs) that the gate rejects on purpose.
- 3.57 M distinct names (2.04 M symbols, 1.52 M alias-only); 81 M indexed named genes; about 22% of all genes have a searchable name.
- Transcripts per gene 1.23, exons per transcript 4.9. A full transcript table for the indexed genes would be about 6 to 8 GB, and every gene about 28 to 34 GB.
- Reading a whole GFF takes 0.1 to 1.8 s, so gene structure can be read on demand for the rows on screen.

## Coverage is the first lever

The cheapest way to make the API more useful is to index more annotations that have real names, not to add namespaces:

- Community annotations with projected gene names (for example TOGA2) are exactly the symbol-searchable, non-model annotations users want. They are skipped today because of a parsing problem (see `TODO-remaining.md`).
- Ab initio GenBank annotations will not match by name. In V1 they are shown as "not searchable by name", not rescued.

## What "not found" means (UI contract)

1. The annotation has no searchable names at all.
2. The annotation is searchable but has no gene with that name.
3. The gene exists but only under a locus-style name that is not indexed.

The UI and the API responses must not collapse these into a single "absent".
