# Project rules

## Memory budget
- Never load large data into RAM. Total use for any script or service step must stay **under about 2 GB**, summed over all parallel workers.
- If a design needs more, choose another strategy instead: streaming with cursors and generators, SQLite scratch tables, external merge sort, smaller batches, fewer workers.
- No `fetchall()` or full in-memory sets/lists over shard-sized or corpus-sized data. Verification code follows the same rule.
- Scratch and temp files go on the NFS work area, never on local disk.
- Measurement scripts report peak RSS; document the expected peak of any new build step.

Why: the host has 11.4 GB of RAM shared with mongod, celery and gunicorn, and swap has been full. A sidecar run that loaded whole shards was killed for lack of memory.
