# Remaining work on the xrefs API refactor

Status as of 2026-11-18. Everything in `API-DESIGN.md` is implemented except the items below.
All code is on main; none of it is deployed until the API image is rebuilt.

## Code, no decision needed
1. **ETag and cache headers.** `ETag` from `data_version` on `/names`, `/hits`, `/loci` (JSON); answer `If-None-Match` with 304. Canonical nginx cache key is a config change (infra side).
2. **Typed `transcript` id in `/sequences`.** Today a transcript must be pinned through a `loci` entry.
3. **Legacy cleanup.** Remove the non-presence paths (`xref_file` `/hits` and `/names`, backend switches) and the `items` form of `/sequences` (kept "one release").
4. **Query replay.** M6 never ran a formal old-vs-new comparison; replay a few hundred real queries.

## Design or decision first
5. **D6, sync.** New, replaced and removed annotations are only picked up by a full rebuild. Also: index `gff_path` / `download_url` go stale after a taxid-folder move until a sync. Largest remaining piece.
6. **D7, namespaces.** Pfam, GO and other non-symbol namespaces (and the annotation filter) are out of scope for v1.

## Needs the maintainer, not code
7. **Deploy.** Rebuild and redeploy the API image. Not live yet: NCBI error classification, missing-GFF lookup across retired-taxid folders, `/hits` (`sources`, `min_matched`, `shape`, TSV, `facets`), prefix terms, `/names` `n_loci` and `sources`, strict parameters.
8. **Annotrieve relocation.** Dry-run, then run, the assembly-file relocation and backfill job on the VM (annotrieve repo changes are still uncommitted).
9. **UI.** See below.
10. **Memory during exports.** Measure peak RSS of a large FASTA export (`docker stats annotrieve-xrefs-api`).

## Coverage: TOGA2 (community registry), investigated 2026-11-19
Sample: Annotrieve annotation `0f20e6413677ad016856197f64bd0ac4` (Cavia aperea, GCA_000688575.1, CommunityRegistry, 13 MB). About 2,040 TOGA2 annotations sit in the registry, roughly 12% of the 17,199 in the corpus.

What the file looks like (why the current harvest drops it):
- **No `Name` attribute at all**, so nothing passes the symbol gate and the annotation is nameless.
- **Comma-joined gene IDs.** `gene` rows carry `ID=hg38.ENSG00000176695,hg38.ENSG00000186092`: one TOGA2 gene is a merged locus of 1 to 4 reference genes (46% have more than one). The `mRNA` rows repeat the same string as `Parent=`. `gene_corpus/stream/attrs.py` splits attribute values on commas, so `ID` becomes two values and children point at genes that do not exist.
- **The symbol lives in the transcript ID.** `ID=<ref>.<ref transcript>#<SYMBOL>#<n>[#paralog]`, for example `hg38.ENST00000585993.3#OR4F17#2511`, `mm10.ENSMUST00000214816#Or4g7#1840`, `HLbosTau10.XM_002690684#OR4K56#823`, `HLbosTau10.XM_025001366#LOC788064#599027#paralog`. The reference prefix varies (hg38, mm10, `HL…` assemblies) and so does the transcript ID style (ENST, ENSMUST, XM_/NM_). One merged gene can have transcripts with different symbols.
- **Row order is not parent-first.** The file is coordinate-sorted, and at the same start `CDS` and `exon` rows come before `gene` and `mRNA`.
- **Size and shape.** 25,546 genes, 82,293 mRNAs (3.2 per gene), 859 k exons; 15,933 distinct symbols; about 45 symbols look like IDs (`LOC…`, `ENS…`) and would be rejected by the gate.

Import sketch (a `toga2` profile, parent-string aware):
1. Treat `ID` and `Parent` as opaque strings (no comma splitting) for this profile; group children by the exact `Parent` string, in a per-annotation dictionary (tens of thousands of entries, small).
2. Name a gene by the distinct symbols found in its transcript IDs (field between the first and second `#`, after the reference prefix); apply the existing symbol gate; the `#paralog` suffix marks a projected paralog (candidate for `via=1`, alias-level).
3. Emit each `ENSG…` in the ID as an `ensembl` xref and keep the reference prefix (hg38, mm10, HL…) as provenance: names are projected, not curated for this species.
4. Test `/sequences` on a TOGA2 slice: `gffread` was not available where this was examined, so whether it copes with the comma Parent strings is unverified (the exon and CDS rows use the mRNA ID, which has no comma).
5. Add a small fixture (a few hundred lines of this file) and a test for the profile.

Open question: projected names are good for cross-species symbol queries, but the UI should label them as "projected by TOGA2 from <reference>".

## UI (next focus)
- Unknown parameters now return `422 unknown_parameter`; check every call the UI makes.
- `/hits` default shape stays `compact`; `records`, `format=tsv`, `facets`, `min_matched`, `sources` are available.
- `*` in names is a prefix term; the response `terms` echoes the expansion (`kind`, `names`, `n_names`).
- `/names` returns `n_loci` and takes `sources`.
- `/loci?format=tsv|ndjson` and `/sequence-estimates` / `-plans` / `-jobs` are ready for the export and FASTA flows.
- See `API-DESIGN.md` section 7 for the intended UI changes.
