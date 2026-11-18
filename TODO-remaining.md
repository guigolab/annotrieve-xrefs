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

## UI (next focus)
- Unknown parameters now return `422 unknown_parameter`; check every call the UI makes.
- `/hits` default shape stays `compact`; `records`, `format=tsv`, `facets`, `min_matched`, `sources` are available.
- `*` in names is a prefix term; the response `terms` echoes the expansion (`kind`, `names`, `n_names`).
- `/names` returns `n_loci` and takes `sources`.
- `/loci?format=tsv|bed|ndjson` and `/sequence-estimates` / `-plans` / `-jobs` are ready for the export and FASTA flows.
- See `API-DESIGN.md` section 7 for the intended UI changes.
