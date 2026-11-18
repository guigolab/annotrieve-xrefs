# Xrefs corpus: speed-ups and new tables (brainstorm)

> **Status: historical brainstorm.** The concrete design that came out of it, with the Phase 0 measurements, is in `API-DESIGN.md` (endpoints) and `PLAN-serve-index-and-cleanup.md` (files, build, drop list). Where they differ, those two win. Main updates since this draft:
> - Measured: the DB is 82 GiB on NFS, cold random reads cost about 19 ms, and the compact name index is about 230 MiB.
> - Dropped from this draft: DFS ordinals (use `annotation_key`), a local-disk copy (not available), delta files (the derived build takes minutes).
> - Added: a second `loci` file for gene coordinates, and prefix-expansion-aware `name_dict`.
> - The tables in section 3 below (`annotation_stats`, `transcript`, `genome_seqid`, `name_group`, `gene_order`, `taxon_stats`) are **deferred**; none is needed for the three v1 endpoints.

Scope: `annotrieve-xrefs` API and services over the ~19k-GFF corpus. Constraint: stay within the current RAM cap (now stated as about 2-3 GB; the corpus is on NFS).

This is the original brainstorm. The concrete build and cleanup plan that came out of it is in `PLAN-serve-index-and-cleanup.md`. Everything here was derived from reading the schema, the `/hits` and `/names` code, the nginx config and `REVISED_PLAN.md`. **Nothing was measured**: the corpus DB was not available. Sizes and latencies are guesses to confirm.

## 1. Where the time and RAM go now

1. **Presence lookups are row-based.**
   - `xref_file` stores `(namespace TEXT, accession TEXT, annotation_key)` per row. The namespace string repeats on every row and the accession is stored in full. The file is large, so little of it stays in page cache.
   - `/hits` does a point lookup per name, then filters by taxon in Python: one `annotation_lineage` query per ~900 keys, in two passes (include, then exclude).
2. **`/names` runs `COUNT(DISTINCT annotation_key) ... GROUP BY` over a prefix range of `xref_file`.** For an unscoped prefix the answer is already in `xref_meta` (`n_annotations`, with an `(namespace, n_annotations DESC)` index). This is probably the slowest cheap endpoint.
3. **Each connection has a 256 MiB page cache, and each thread opens its own connection.**
   - `READ_CACHE_SIZE=-262144` is set per connection in `session.py`. Connections are thread-local and sync FastAPI handlers run in a threadpool of about 40 threads per worker. With 2 gunicorn workers the worst case is about 80 x 256 MiB of private cache holding the same pages.
   - `mmap_size=0` is set on purpose, so the OS page cache cannot be shared between threads and workers.
   - Check this first against the RAM cap.
4. **No precomputation for rollups.** Any "how many annotations per taxon" question is a join against `annotation_lineage`.

## 2. Speed-ups that need no extra RAM

- **Shrink presence into per-name postings.** With ~19k annotations, a name's presence is a sorted list of annotation ids in a small blob.
  - One table: `name_posting(name_id INTEGER PK, n INTEGER, postings BLOB)`, delta-varint encoded. A name in a few annotations costs a few bytes. A very common name costs about 2.4 KB as a bitset.
  - A `name_dict(name_id, name_cf)` for the names.
  - Guess: 5-10x smaller than `xref_file` for symbol and alias. Needs a measurement.
- **Make the taxon filter a bit operation.**
  - Assign annotation ordinals in taxonomy DFS order (as `REVISED_PLAN.md` proposes for `genome_ord`). A clade is then a contiguous ordinal range.
  - Store `taxon_range(taxid, lo, hi)`, which is tiny, and optionally a bitmap for any non-contiguous scope.
  - With Python ints as bitsets, a search is: fetch N blobs, OR/AND them, mask by the include range, clear the exclude ranges. That is milliseconds, with no lineage join and no 900-variable chunking.
  - Distribution rollups per child taxon become popcounts over ranges (about 100 children x 200 names is on the order of 20k small operations). No precomputed taxon x name table is needed.
- **Serve `/names` from precomputed counts.** Unscoped: from `xref_meta`. Scoped: range-scan the name dictionary for candidates, then popcount each posting under the scope range.
- **Use `json_each(?)` or a VALUES CTE** instead of chunked `IN (?,?,...)` (already in the plan).
- **Cap and share memory:**
  - Limit the threadpool (about 8 threads).
  - Set `cache_size` to about 16-32 MiB per connection.
  - Turn on `mmap_size` on the immutable read-only DB so workers and threads share one copy through the OS page cache, which the kernel can reclaim. **On NFS, only do this once the file is on local disk** (see section 4).
- **Build-time layout:** insert in primary-key order, 16-32 KB `page_size`, one `ANALYZE`, covering indexes, no TEXT in hot keys.
- **Cache better:**
  - nginx caches GET for 24 h, but the key is `request_uri`. `names=BRCA1,TP53` and `names=TP53,brca1` are different entries. Normalize (lowercase, sort, dedupe) on the client, or on the server with a canonical redirect.
  - POST is not cached at all. Add a small byte-bounded in-process LRU keyed by a hash of the normalized request plus a `data_version`.
  - Run a warm-up script after each corpus publish (example gene sets plus top queries from nginx logs).

## 3. New tables worth building

All can be computed in a pass over the corpus.

| Table | What it gives you |
|---|---|
| `name_dict` + `name_posting` | The fast presence path above |
| `name_norm` column in `name_dict` (casefold, strip `-_ .`, normalize Greek letters) | Exact lookups that survive `IL-6` vs `IL6` vs `Il6`, and a cheap "did you mean" with no FTS |
| `annotation_stats(annotation, n_genes, n_coding, n_named, named_fraction, n_transcripts, mean_exons, busco)` | The honesty layer: "this annotation names only 12% of its genes, so a missing symbol means little." Also feeds the representative-genome rule in the plan and a quality filter |
| `transcript(genome, local_id, transcript_id, start, end, n_exons, cds_len, is_longest)` | `transcript: "longest"` resolved in SQL, FASTA size estimates, isoform lists without parsing a GFF. The biggest win for the FASTA service |
| `genome_seqid(..., ncbi_accession, length)` | Resolves seqid to accession and length without reading assembly reports at request time. Lets `/estimate` know the bp to fetch |
| `gene.voffset` (bgzf virtual offset), already in the plan | Slice a gene's GFF lines directly instead of a tabix lookup per gene |
| `name_group(name_id, group_id, n)`, built from genes where a symbol and an `ncbiortholog` (or `geneid`) co-occur | Synonym expansion across species, which addresses the biggest scientific caveat (the same gene has different names in different taxa). The UI can offer "also include these 6 names". `ncbiortholog` is already in the Tier A namespaces |
| `gene_order(genome, seqid_ord, rank, local_id)` | Left and right neighbours of a gene; cheap to add microsynteny later |
| `taxon_stats(taxid, n_annotations, n_organisms)` | Instant denominators for "found in 62% of Mammalia". Tiny |

For FASTA, the most useful are `transcript` and `genome_seqid`: the first removes client-side transcript resolution and enables a real dry-run `/estimate`, the second removes per-request file reads.

## 4. NFS-specific cautions (added after the first draft)

- The corpus is on NFS and each B-tree page miss is a network round-trip, so copy the serving DB to a local volume at publish time. This is probably the biggest latency win.
- Mapped NFS files can raise `SIGBUS` if truncated or replaced underneath. Publishing by rename to a new inode is safer, but enable mmap only after the file is local.
- Per-annotation shard opens (19k files) stay on NFS. A local LRU of hot shards is an optional later step.

## 5. Suggested order

1. **Measure first** (`REVISED_PLAN.md` phase 0 asks for this). On the real DB:
   - `dbstat` sizes for `xref_file`, `annotation_lineage` and `gene_hit`.
   - `EXPLAIN QUERY PLAN` for a 200-name `/hits` and a broad `/names`.
   - Cold vs warm timing.
   - Peak RSS under about 20 concurrent requests.
2. **Fix memory:** thread cap, smaller per-connection cache. Configuration only, no rebuild.
3. **`/names` from `xref_meta`** for unscoped prefixes.
4. **Build `name_dict`, `name_posting` and DFS ordinals**, then switch `/hits` to bitset operations, keeping the old path behind a flag to compare results.
5. **`transcript` and `genome_seqid`**, then the FASTA `/estimate`.
6. **`annotation_stats`, then `name_group`**, since these carry the most scientific value.

## 6. Open questions from the first draft

- Does `REVISED_PLAN.md` (the `genome` / `organism` / `gene.db` design) supersede the live `xref_file` path? This brainstorm assumed the postings idea sits alongside the plan's `genome` table, not instead of it.
- What is the actual RAM cap, and is the corpus on local disk or NFS? (Since answered: about 2-3 GB, and NFS.)
- How big is `xref_file` today? If it is a few GB, it fits in the page cache anyway and the postings change would be mostly about CPU.
