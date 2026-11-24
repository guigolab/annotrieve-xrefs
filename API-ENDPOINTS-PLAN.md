# API endpoints: consistency fixes and the sequences plan

Status: proposal written 2026-11-24, nothing implemented yet. Tick the boxes in section 4 as work is done and add a line to the progress log (section 5).

Scope: the public routes of the xrefs API (`/hits`, `/names`, `/loci`, `/annotations/...`, `/sequences...`). Memory and speed decisions are in `API-LOCKED-DECISIONS.md`; this file is only about how the routes look and behave.

**Compatibility rule for everything below (the sequences change is the exception: nothing uses `/sequences` yet except the genes page, which is still in development, so job mode is removed outright, including the `mode` and `sync_max_genes` estimate fields):** a change may add a name, a field or a header, but never removes or renames what clients use today. Old names stay as aliases and answer with a `Deprecation: true` header, so scripts keep working. The only intentional removal is the job mode of `/sequences` (section 2), which you asked for.

---

## 1. Consistency fixes

Every "today" below was checked on 2026-11-24 in the app's own route list, in the code, or on the live API.

### 1.1 Paging works three different ways

| Route | Today |
|---|---|
| `GET /hits` | no paging: always returns every match (up to about 18,000 rows, 9 MB) |
| `GET /loci` | `limit` + `cursor`, answer has `next_cursor` |
| `GET /annotations/{id}/genes`, `.../accessions`, `.../xrefs` | `limit` + `offset`, answer has `total`, `limit`, `offset` |
| `GET /annotations/{id}/namespaces` | answer has `total`, `limit`, `offset` although it never pages (`limit` simply equals `total`) |

Problem: a client needs three loops. The `namespaces` envelope promises paging that does not exist.

Fix:
- Every route that pages accepts `limit` and `cursor`, and answers with `next_cursor` (null on the last page). `/annotations/*` keeps `offset` as an alias and `total`/`offset` in the answer, and adds `next_cursor`. The cursor is the encoded offset there, so no query changes.
- `/hits` stays unpaged on purpose: it is one complete answer (the `complete` guarantee, the TSV export and the facets all depend on that). This is written into the route description with the real maximum size, so nobody has to guess.
- `/annotations/{id}/namespaces` keeps its fields for old clients but they are marked deprecated in the route description; new clients read `results` only.

### 1.2 A gene has two ids

`/annotations/{id}/genes/{local_id}` takes `local_id`. `/loci` returns, and `/sequences` accepts, `locus_ref` (`annotation_id:local_id`).

Fix: the annotation is already in the nested path, so the nested routes keep `local_id`. Gene cards and `/loci` rows keep `locus_ref`, and a new route `GET /genes/{locus_ref}` (and `/genes/{locus_ref}/xrefs`) answers the same as the nested ones. The nested routes stay. One id then works everywhere, including in a URL you copy from a `/loci` answer.

### 1.3 Similar filters have different names or are missing

| Today | Fix |
|---|---|
| `sources` exists on `/hits` and `/names`, not on `/loci` | add `sources` to `/loci` (same values: Ensembl, GenBank, RefSeq) |
| `/names` takes `prefix`; the accessions route takes `q` for an accession prefix | the accessions route also accepts `prefix` (`q` stays as alias) |
| `/annotations/{id}/genes` has `names`, `curies` and `match` | unchanged (they mean different things: name list, accession list, and any/all); each is described in plain words in the route text |

### 1.4 Two output switches: `shape` and `format`

`shape` is `compact` or `records`; `format` is `json` or `tsv`. TSV silently means records, whatever `shape` says, and `shape` shows a default of `compact` in the route list.

Fix: `shape` has no default value in the route list (null means "the usual one for this format"). `format=tsv` with an explicit `shape=compact` answers 422 `invalid_combination` instead of quietly ignoring it. The answer echoes `shape` so the client sees what it got.

### 1.5 Two status codes for bad input

Validation problems answer 422, other query errors answer 400, and the table that decides which is copied into three route files (`hits`, `names`, `annotations`) with a fourth variant in `loci`.

Fix: one function in one file decides; every input problem answers 422, with no exception: a bad cursor was already 422 in `/loci`, so 400 is not used at all any more. The copies are deleted. A test lists every error code the API can raise and checks it has exactly one status.

### 1.6 Different key names for lists

`/hits` compact uses `rows` (positional, described by `fields`); `/hits` records, `/names`, `/loci` and `/annotations/*` use `results` (objects). That is a real difference, not an accident: positional rows are about half the size.

Fix: no rename. It is written down in the route text: `rows` means "list of arrays, see `fields`", `results` means "list of objects". The tests check every list answer has one of the two, never both.

### 1.7 Namespaces

Today `/annotations/{id}/namespaces` returns 13 entries for the human RefSeq annotation. Eleven are cross-reference databases (dbsnp, genbank, geneid, hgnc, imgt/gene-db, mim, mirbase, refseq_protein, refseq_transcript, rfam, vista). Two are not databases: `symbol` and `alias`, the name types `/hits` searches. A name like `imgt/gene-db` has a slash, which forces the odd `{namespace:path}` route.

Fix:
- Each entry gets `kind` (`database` or `name`) and `slug` (`imgt-gene-db`: lowercase, slashes become dashes). The accessions route accepts the slug, the old literal name, and the percent-encoded form.
- The route text says what namespaces are for: finding genes by an accession inside one annotation. There is no search by accession across annotations; that is stated, not hidden. (Adding one would be a separate, larger feature.)
- `symbol` and `alias` stay in the list for old clients but clients can filter on `kind`.

### 1.8 Things I considered and would leave alone

- Route nouns (`/hits`, `/names`, `/loci`): they are search endpoints, not resources, and renaming them breaks every client for no gain.
- POST for long name lists on `/hits` and `/loci`: a normal, accepted pattern.
- The summary table for `/annotations/{id}/namespaces` (speed only): dropped at your request.

---

## 2. Sequences: synchronous only, at most 200 sequences

### Today (verified)

- `POST /sequences` has two ways to answer. Up to `SEQUENCES_SYNC_MAX_GENES` = 50 genes the FASTA comes back in the response. Above that (up to `SEQUENCES_MAX_GENES` = 2,000) it answers 202 and starts a background job that the client polls at `GET /sequence-jobs/{id}`, then downloads from `GET /sequences/results/{id}`.
- Finished results are kept on the network drive (24 h, at most 100) so an identical request is served again at once, and the per-gene outcome table lives at `.../manifest.tsv`.
- Per client IP: at most 2 exports running and 20 per hour.
- nginx waits up to 900 s for these routes.
- Live timing on 2026-11-24 (Ensembl human annotation, `transcript`): 10 genes took 3.0 s, 50 genes took 13.4 s, so about 0.27 s per gene, almost all of it fetching sequence from NCBI. A request of 200 would therefore take about 55 s. This is an extrapolation from two points; it is measured again after the change (section 4).
- Found while measuring: on the live API the human RefSeq annotation `a63cdc06...` answered every gene with `skipped: gff_unavailable` ("annotation GFF file not found on disk"), while the Ensembl and GenBank human annotations worked. That is a data or mount problem on the server, not a code problem, and it is listed in section 4 as something to check there.

### Decision

One mode: the FASTA is in the response, or the request is refused. At most **200 sequences** per request.

| Part | Change |
|---|---|
| Limit | new setting `SEQUENCES_MAX_SEQUENCES` = 200. Genes are counted first from the index (cheap); more than 200 answers 413 `too_many_genes` as today. `max_genes` in the body can only lower it. |
| Sequences versus genes | one gene gives one sequence, except with `transcript: "all"` or a rule that keeps several transcripts. So the genes are resolved first (local GFF reads, no network), the sequences are counted, and a total above 200 answers 413 `too_many_sequences` with the count, **before** any NCBI call. |
| Job mode | removed: the 202 answer, `Location` header, `GET /sequence-jobs/{id}`, the `SequenceJob` model, the background thread code, `SEQUENCES_SYNC_MAX_GENES`, `SEQUENCES_MAX_GENES` (2,000). |
| Kept | `POST /sequence-estimates` and `POST /sequence-plans` (both cheap and local; the plan already caps at 200 genes, the same number), the stored results with `.../results/{id}` and `.../manifest.tsv` (the response points to the manifest, and repeats are free), the per-IP limits, the memory limits from the other document. |
| Time limit | `SEQUENCES_DEADLINE_S` 800 to 120 (about twice the expected 55 s). When it runs out: 504 `sequences_timeout`, no partial FASTA, the client narrows the request. |
| Concurrency | stays at 1 export per worker, so 2 at the same time in the whole service. A third gets 503 `sequences_busy` with `Retry-After`. With 55 s per full export that is plenty for about 300 users a month; it is revisited if the logs show busy answers. |
| nginx | read timeout for `/sequences` 900 s to 180 s, so a stuck connection is freed sooner. |
| Estimate answer | `mode` and `sync_max_genes` go away (there is only one mode); `max_sequences` is added. Old clients reading `mode` get the constant `"sync"` for one release, then it is dropped. |

### What users see

| Request | Answer |
|---|---|
| up to 200 sequences | `200` with the FASTA; headers `X-Sequence-Count`, `X-Skipped-Count`, `X-Manifest-Url`, as today |
| more than 200 genes or sequences | `413` with the count found and the limit |
| too many exports now | `503 sequences_busy` or `429` (per IP), with `Retry-After` |
| took longer than 120 s | `504 sequences_timeout` |

### Risks

- A browser or script that waits less than about a minute will time out on a 200-sequence request. The estimate route gives the gene count and `est_bp_upper_bound` first, and the route text says to use it for big requests.
- Any client that handles the 202 answer (the web front end, `annocli`) must be searched before the code is removed. That search is the first step in section 4.

---

## 3. Order of work

1. Sequences (section 2): it is the only change that removes something, and it makes `/sequences` simpler.
2. 1.5 status codes (one decision function), then 1.4 `shape` and `format`.
3. 1.1 paging and 1.3 filters.
4. 1.7 namespaces `kind` and `slug`, then 1.2 the `/genes/{locus_ref}` routes.

Each step ends with: tests, the committed `api/openapi.json` regenerated, a line in the progress log, commit and push to `main`.

## 4. Work list

- [x] Find every client of `/sequence-jobs` and of the 202 answer: only the genes page of `annotrieve-genes-ui` (`front/lib/api/xrefs.ts`, `xrefs-types.ts`, `xrefs-sequences.test.ts`); nothing else uses sequences yet (confirmed by the owner), so no compatibility is kept
- [x] Sequences: 200-sequence limit, resolve before fetch, remove job mode, `SEQUENCES_DEADLINE_S` 120, estimate fields, nginx timeout
- [ ] Sequences: measure a full 200-sequence request on the real stack and write the time here
- [ ] Check on the server why the human RefSeq annotation `a63cdc06...` answers `gff_unavailable`
- [x] 1.5 one status-code decision (`helpers/errors.http_for_query_error`, always 422; a test fails if any route file contains `status_code=400`)
- [x] 1.4 `shape` without default, echoed in the answer, 422 `invalid_combination` on tsv + compact
- [x] 1.1 `cursor` on `/annotations/*` (genes, accessions, xrefs), route text for `/hits` and `/namespaces`
- [x] 1.3 `sources` on `/loci` (GET and POST, echoed as `source_filter`, part of the cursor key), `prefix` on the accessions route
- [ ] 1.7 `kind` and `slug` for namespaces
- [ ] 1.2 `/genes/{locus_ref}` routes
- [x] 1.6 route text for `rows` and `results`, test that every list has exactly one
- [x] Regenerate `api/openapi.json` (the committed one was stale: its test failed before this work; regenerated inside the production image)

## 5. Progress log

- 2026-11-24: plan written; live sequences timing taken (10 genes 3.0 s, 50 genes 13.4 s); the `gff_unavailable` finding recorded.
- 2026-11-24: sequences done. One mode only: `200` with the FASTA. Limits `SEQUENCES_MAX_GENES` 200 and `SEQUENCES_MAX_SEQUENCES` 200 (the gene limit is also capped by the sequence limit); when the transcript rule can keep several per gene the sequences are counted from the GFF slices before any NCBI call (`413 too_many_sequences`), with a second guard inside the export; identical request in progress: `409 export_running`; deadline 120 s; nginx read timeout 180 s for `/sequences` and `/sequence-plans`. Removed: 202, `Location`, `GET /sequence-jobs/{id}`, `SequenceJob`, `start_job` and its thread, `SEQUENCES_SYNC_MAX_GENES`, estimate fields `mode` and `sync_max_genes` (added `max_sequences`). Kept: estimates, plans, stored results with manifest, per-IP limits. Tests: 379 api (including the openapi test, fixed by regenerating the file) and 221 gene_corpus pass. Not measured yet: a full 200-sequence request on the real stack (needs NCBI and the real GFF files). The genes page still handles a `202` and polls `/sequence-jobs`; that code in `annotrieve-genes-ui` is now dead and can be deleted.
- 2026-11-24: 1.5 and 1.4 done. All four copies of the code table are gone; `invalid_offset` and `prefix_too_broad`-style errors that answered 400 now answer 422. `/hits` answers carry `shape`; `shape` has no default in the route list; `format=tsv` with `shape=compact` answers 422. `api/openapi.json` regenerated in the production image. Tests: 384 pass.
- 2026-11-24: 1.1, 1.3 and 1.6 done. The `/annotations/*` lists take `cursor` and answer `next_cursor` (null on the last page and on `/namespaces`); the cursor is the encoded offset, so the depth limit (10,000 rows) is unchanged and is stated in the error; giving both `cursor` and `offset` answers 422 `invalid_combination`. `/loci` takes `sources` like `/hits` and `/names`. The accessions route takes `prefix` (same as `q`; different values: 422). Route text for `/hits` (never paged, size, `rows` versus `results`) and `/namespaces` (never paged, what `symbol` and `alias` are). `api/openapi.json` regenerated. Tests: 391 pass. Left: 1.7 namespaces `kind` and `slug`, 1.2 `/genes/{locus_ref}`.
