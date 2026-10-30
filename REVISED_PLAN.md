# Annotrieve gene-set search: implementation plan

Design decisions and build order for the xref sidecar, as agreed on 18 Oct 2026.

## 1. Goal and locked decisions

Users search a gene set (names, aliases, prefixes) under a taxon (Eukaryota by default) and want either a list of annotations or gene rows. The sidecar complements Annotrieve rather than replicating NCBI/Ensembl.

- **Term vs key.** A *term* is what the user typed (`seleno*`). A *key* is a `(namespace, accession)` the term expands to. A genome matches a term if it has **any** of its keys. Coverage is counted over terms.
- **Unit of answer is the genome**, meaning one row per (canonical assembly, source) using the GCF↔GCA pairing. Genomes roll up into an organism, which is the row users see first (section 2.3). Individual assemblies and annotations are a drill-down.
- **Coverage, not just any/all.** Every row reports `n_matched / n_terms` and the missing terms.
- **Zoom out before paginating.** The default result is a rollup by child taxon, and genomes are listed only once the scope is small or the user asks.
- **Species and clades are the same thing:** a taxon node, scoped through `annotation_lineage`.
- **Prefixes** are resolved before the query. The user can keep a prefix as one OR term, untick keys, or split it into one term per key. Counts for a prefix are unions, never sums.
- **Gene rows are key-major** and come from a dedicated narrow-plus-wide pair of tables, in a separate database file.
- **Per-annotation GFF/csi files stay the source of truth**, used only for subset-GFF extraction.

**Organisms and multiple assemblies.** A species can have hundreds of assemblies (Homo sapiens: 474 assemblies, 535 annotations). The default row is therefore an **organism**, scored by its representative genome: coverage is the best single genome, never a union across assemblies. The row also shows how many assemblies hold each gene (for example 470/474), and assemblies with identical results are grouped instead of listed one by one.

## 2. Data model

### 2.1 `xref.db` (presence and coverage)

Keep `xref_meta` and `xref_file` as they are, and add:

```sql
CREATE TABLE genome (
    genome_ord     INTEGER PRIMARY KEY,   -- dense; assigned in taxonomy DFS order
    assembly_canon TEXT NOT NULL,         -- GCF/GCA paired canonical accession
    source         INTEGER NOT NULL,      -- 1 Ensembl, 2 GenBank, 3 RefSeq, 4 community
    taxid          INTEGER NOT NULL,
    annotation_key INTEGER NOT NULL,      -- current representative (latest)
    n_versions     INTEGER NOT NULL,
    UNIQUE (assembly_canon, source)
);
```

Also add to `xref_meta`: `n_genomes` (distinct genomes with the key), next to the existing `n_annotations`. Add n\_organisms as well; the UI shows organisms first, then genomes.

Placeholder rule at import: names matching `LOC\d+`, `gene-*`, and unnamed Ensembl genes are **not indexed as keys**. Decide the exact patterns once and keep them in one function.

### 2.2 `gene.db` (gene rows, separate file, attached read-only for exports)

```sql
CREATE TABLE key_dict (
    key_id INTEGER PRIMARY KEY,
    namespace TEXT NOT NULL,
    accession TEXT NOT NULL,
    UNIQUE (namespace, accession)
);

CREATE TABLE gene_hit (                    -- narrow, about 10 bytes per row
    key_id     INTEGER NOT NULL,
    genome_ord INTEGER NOT NULL,
    local_id   INTEGER NOT NULL,
    PRIMARY KEY (key_id, genome_ord, local_id)
) WITHOUT ROWID;

CREATE TABLE gene (                        -- one row per gene, not per key
    genome_ord   INTEGER NOT NULL,
    local_id     INTEGER NOT NULL,
    seqid_ord    INTEGER NOT NULL,
    start        INTEGER NOT NULL,
    end          INTEGER NOT NULL,
    strand       INTEGER NOT NULL,
    feature_type INTEGER NOT NULL,
    biotype      INTEGER,
    primary_name TEXT,
    voffset      INTEGER,                  -- bgzf virtual offset, for subset GFF (phase 4)
    PRIMARY KEY (genome_ord, local_id)
);

CREATE TABLE genome_seqid (genome_ord INTEGER, seqid_ord INTEGER, seqid TEXT,
                           PRIMARY KEY (genome_ord, seqid_ord)) WITHOUT ROWID;
CREATE TABLE dict_feature_type (id INTEGER PRIMARY KEY, name TEXT UNIQUE);
CREATE TABLE dict_biotype (id INTEGER PRIMARY KEY, name TEXT UNIQUE);
```

Scope rules: representative annotation per genome only, gene features only, named keys only. Older versions are not in `gene.db`; for those the user gets the annotation ID and goes to Annotrieve.

Pragmas on read connections: `page_size` 16–32 KB (set at build time), large `cache_size`, `mmap_size`, WAL, `query_only`.

### 2.3 Organism and representative genome

```sql
CREATE TABLE organism (
    organism_ord   INTEGER PRIMARY KEY,
    taxid          INTEGER NOT NULL UNIQUE,   -- definition pending, see section 8
    n_genomes      INTEGER NOT NULL,
    rep_genome_ord INTEGER NOT NULL
);
-- genome gains: organism_ord, is_reference, busco_complete, rep_reason
```

Assign `genome_ord` in taxonomy DFS order and then by organism, so an organism is a contiguous `genome_ord` range.

The representative genome is picked at import with one rule, and the reason is stored in `rep_reason` so users can see why it was chosen:

1. Reference or default assembly (RefSeq reference, Ensembl default).
2. Highest BUSCO complete, with a missing value counting as lowest.
3. Latest release. Remaining ties go to RefSeq, then Ensembl, then GenBank.

Gene-count sanity filters (for example partial annotations with a very low non-coding count) are not part of v1. `gene.db` still stores the latest annotation of every genome, so an all-assemblies gene export stays possible; if its size becomes a problem, store representatives only.

## 3. Query core

Every endpoint shares one pipeline. Resolved keys are passed per request (`json_each(:keys)` or a `VALUES` CTE), never through shared temp tables.

1. **Resolve:** terms → keys, with `term_idx` per key (after the user's untick and split choices).
2. **Genome-first collapse** (one row per genome and term, which dedupes versions and multi-key hits):

   ```sql
   SELECT a.assembly_canon, a.source, q.term_idx
   FROM q
   JOIN xref_file xf ON xf.namespace=q.namespace AND xf.accession=q.accession
   JOIN annotation a ON a.annotation_key = xf.annotation_key
   WHERE a.is_latest = 1
   GROUP BY a.assembly_canon, a.source, q.term_idx;
   ```
3. **Coverage per genome:** `COUNT(*)` of that result grouped by genome, plus the matched term set (`GROUP_CONCAT(term_idx)`).
4. **Scope and rollup:** only now join `annotation_lineage`. Scope is `include_taxids` / `exclude_taxids`; for `group_by=child`, map each genome to the child of the scope that contains it. Eukaryota with no scope skips the scope filter and uses the join only to find children.
5. **Best source per assembly:** take the max `n_matched` across sources when a row per assembly is requested.

6) **Organism step:** aggregate the per-genome coverage to organisms.

```sql
SELECT g.organism_ord,
       COUNT(*)                    AS n_genomes_hit,
       MAX(c.n_matched)            AS best_cov,
       SUM(c.n_matched = :n_terms) AS n_full
FROM coverage c JOIN genome g USING (assembly_canon, source)
GROUP BY g.organism_ord;
```

Per-term support is a `COUNT(DISTINCT assembly_canon)` per (organism, term) over the same genome-first result. Signature groups (source plus matched terms) inside one organism are a `GROUP BY` on the per-genome matched set. Taxon rollups count organisms by default, with a switch to count genomes.

The cost is proportional to matched rows (roughly 80 keys × 4k annotations ≈ 320k), not to the 475M-row table.

## 4. API

### `POST /resolve`

Input: raw terms plus scope. Output, per term: candidate keys with `n_genomes` **in scope**, union count for the term, a truncated flag, and an ambiguity or zero-hit flag. Prefix expansion is a range scan (`>= 'seleno' AND < 'selenp'`, never `LIKE`), with a minimum prefix length of 3 and a cap of about 500 keys that also returns the true count.

### `POST /hits` (`GET` keeps the old shape as an alias)

Body: `scope`, `terms` (resolved keys, `split` flags), `filters` (`sources`, `latest_only`, `min_coverage`, `match`), `group_by` = `child | organism | genome | annotation`, plus pagination.

- `child`: one row per child taxon with `n_genomes`, `full`, `mean_cov`, per-term genome counts and a missing-frequent list.
- `genome`: one row per assembly and source with the representative annotation, `n_matched`, matched and missing terms, and one cell per source. Identical result rows may be collapsed with a count.
- `annotation`: every version for a chosen genome (the expand action).
- Return facets (by source, by rank) with every response.

`organism` is the new default list view: one row per organism with its representative genome, `best_cov`, `n_genomes`, `n_full` and per-term support (`k/n` assemblies). Expanding an organism returns signature groups (assemblies with the same source and matched terms, a count and one example), then individual assemblies, paginated. `child` rollups count organisms; pass `unit=genome` to count genomes instead.

For `/hits`, send an `estimate` before heavy queries and refuse or warn when a broad prefix would match too many rows (use `xref_meta` counts).

### `POST /exports`

```json
{"kind": "matrix|genes|annotations", "scope": {...}, "terms": [{"label": "seleno*", "keys": [...], "split": false}],
 "filters": {...}, "grain": "genome|annotation",
 "selection": {"exclude_assemblies": [...]}, "format": "tsv|ndjson|bed",
 "layout": "wide|long"}
```

Send resolved keys and a filter plus an exclude list, never thousands of selected IDs.

`grain` is `organism` (default), `genome` or `annotation`.

## 5. Exports

### Delivery

- The request computes an estimate (matrix: genomes × terms; genes: matched genome-key pairs as a lower bound).
- **Small** (about < 200k cells or < 50k gene rows): stream with `Content-Disposition`, so `curl` works.
- **Large:** `202 {job_id, estimate}`; `GET /exports/{id}` for status and progress; `GET /exports/{id}/download`. Over the hard cap: `413` with the estimate and a hint to narrow the scope.
- **Cache** by hash of the normalized request plus the data version (last import id). Files expire after a TTL (24 h). Write to a temp file, gzip, then atomically rename.
- **Limits:** N concurrent jobs and one per client. Every file starts with `#` provenance lines: query, scope, key expansions, data version, timestamp.

### Matrix (from `xref.db` only)

Step 2 of the query core, then `GROUP_CONCAT` the present terms per genome and expand into a row in the stream writer, so memory stays constant. Wide layout by default and long layout for R/pandas. Columns sorted by user order or prevalence, with a cap of about 2,000 terms. Cell modes: `presence` now; `count` only once gene-level data is available (see phase 3).

### Annotations

A list of annotation IDs with metadata (taxid, species, assembly, source, version, `n_matched`, missing terms), ready to feed into Annotrieve.

### Genes (from `gene.db`)

1. Take the selected genomes (the `sel` set) from the genome-first step.
2. Per key, range-scan `gene_hit` and keep rows whose `genome_ord` is in `sel`. If `genome_ord` is in taxonomy DFS order, clade scope is a contiguous range.
3. Look up `gene` by `(genome_ord, local_id)` in batches sorted by that pair.
4. Deduplicate on `(genome, gene, term)`, since one gene can match via symbol and alias.
5. Output: `taxid, species, assembly, source, annotation_id, term, matched_key, namespace, gene_id, symbol, seqid, start, end, strand, biotype`. BED output converts to 0-based half-open and puts `assembly|gene_id|symbol` in the name column.

Gene rows always come from the same representative annotation as the matrix row.

### Organism grain and multi-assembly species

- **Matrix:** one row per organism using its representative genome, plus a support column per term (`k/n` assemblies). `grain=genome` returns every assembly.
- **Gene rows:** default to the representative genome of each organism. An all-assemblies export is opt-in and shows its estimate first, since human alone adds 474 rows per key.
- **Annotation list:** the representative annotation per organism by default, and every version with `grain=annotation`.

## 6. Import pipeline changes

- Maintain the `genome` table: assign `genome_ord` in taxonomy DFS order, set the representative on each import, and bump `n_versions`.
- Apply the placeholder rule and the gene-feature filter once, and feed both `xref_*` and `gene.db` from the same pass so they cannot drift.
- Build `gene.db` per genome in one transaction. When a new version replaces the representative, delete the old rows by full primary key (the old gene list is known) or rebuild the genome's rows, then insert.
- Keep `gene.db` rebuildable from the GFFs without touching `xref.db`.
- Bump a `data_version` on every import (used for export cache keys).

Compute `organism` and each representative genome in the same import pass (reference flag, BUSCO, release date), and recompute an organism's representative whenever one of its genomes changes.

## 7. Build order and exit criteria

**Phase 0: measure and verify (1–2 days)**

- Confirm `annotation_lineage` stores all ancestors (if leaf-only, build a closure table).
- Confirm `annotation` has canonical assembly, source and an `is_latest` flag; add them if not.
- On the old `gene_hit`: row count, file size, `dbstat`, `EXPLAIN QUERY PLAN` of an 80-key query, cold versus warm timing.
- *Exit:* the genome-first query core runs in acceptable time on current tables for Eukaryota with 3 broad prefix terms. If it doesn't, move the presence tier to `xref_genome(key_id, genome_ord)`, latest-only, or bitmaps.

* Check the genomes-per-organism distribution, and confirm the reference flag and BUSCO exist for most annotations (human: 534 of 535 have BUSCO).
* Decide the organism definition (section 8).

**Phase 1: grouped search on current tables**

- `genome` table, `/resolve` (scoped counts, union counts, caps), `/hits` with `group_by=child|genome`, coverage and missing terms, facets.
- *Exit:* the Eukaryota → clade → genome drill-down works end to end in the UI.

* `organism` table, the representative-genome rule, and `group_by=organism` with per-term support and signature groups.

**Phase 2: matrix and annotation-list exports**

- `/exports` with request shape, estimate, streaming, caching and provenance; wide and long layouts.
- *Exit:* a 5k-genome × 80-term matrix streams or completes as a job within the agreed time.

* `grain=organism` as the default for the matrix and the annotation list.

**Phase 3: gene rows**

- Build `gene.db` (new schema, placeholders skipped, representative only), then `kind=genes` in TSV/NDJSON/BED with job and cache support. Optionally add per-genome copy number so the matrix can offer `count`.
- *Exit:* an 80-key, 1k-genome gene export on a cold cache within the agreed time; measure `gene.db` size against the estimate.

* Gene rows default to the representative genome; all-assemblies is opt-in with an estimate.

**Phase 4: subset GFF3 (only if users ask)**

- Use `voffset` to read gene models from the bgz and write a zip with one file per annotation.

**Later:** cross-source concordance on the same assembly (needs seqid mapping through the assembly report), stable-ID namespaces from Dbxrefs (GeneID, HGNC, Ensembl IDs), ortholog-based synonym expansion, and TOGA2 using the same `genome` logic.

## 8. Open items and risks

- **Row counts of the new `gene_hit`.** The expected size reduction comes from guesses, so verify it in phase 0 with a sample import.
- **Placeholder patterns** need a decision (and a test set), because they drive table size more than anything else.
- **Broad prefixes** (`se*`) need the estimate-and-refuse path in both `/resolve` and `/hits`.
- **Alias collisions across species** mean results are "annotated under this name", not proof of the gene's presence. Report the matched namespace and treat alias hits as lower confidence.
- **Quality metadata** for `representative=best` (BUSCO, gene count, reference flag) is not available yet, so "latest" is the only rule for now.
- **The 500k–1.5M-loci outliers:** check what they are (fragmented assemblies or non-gene features) and cap them at import if needed.

* **Organism definition:** is it the taxid assemblies are assigned to (strains and subspecies can be their own organism), or always a species? If strains roll up into species, derive it through the taxonomy.
* **Representative rule inputs:** confirm the reference flag and BUSCO are stored per annotation. Outlier checks on gene counts are deferred.
* **Coverage across assemblies:** organism coverage is the best single genome. A union across assemblies would only be shown as a flag, if at all.
