"""
Per-annotation sidecar: a small sorted file written from one shard.

Layout (``<per_gff>/<annotation_id>/sidecar.tsv.gz``, gzip text, see
``PLAN-build-tables.md`` section 3)::

    #sidecar_version=1
    #annotation_id=...
    ... more ``#key=value`` header lines ...
    #names
    name<TAB>kind<TAB>local_id          (kind: s = symbol, a = alias)
    #genes
    local_id<TAB>source_gene_id<TAB>seqid<TAB>start<TAB>end<TAB>strand<TAB>feature_type<TAB>biotype<TAB>primary_name
    #end names=<n> genes=<n>

* ``#names`` is sorted by (name, kind, local_id); names are casefolded.
* ``#genes`` is sorted by ``local_id`` and holds **only genes that have at least one
  name row in ``#names``** (unnamed genes cannot be reached by name; the shard keeps
  every gene). NULL text values are written as an empty field.
* The trailing ``#end`` line lets readers detect a truncated file.
* The header stores a fingerprint of the shard it was built from. Shards are
  rewritten in place, so staleness is detected by fingerprint, never by path.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import os
import sqlite3
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SIDECAR_VERSION = 1
SIDECAR_NAME = "sidecar.tsv.gz"
SHARD_NAME = "genes.sqlite"

KIND_SYMBOL = "s"
KIND_ALIAS = "a"
_NAMESPACE_KIND = {"symbol": KIND_SYMBOL, "alias": KIND_ALIAS}

_BAD = str.maketrans({"\t": " ", "\n": " ", "\r": " ", "\0": ""})


class SidecarError(Exception):
    """Sidecar missing, truncated, malformed or of an unsupported version."""


@dataclass(frozen=True)
class SidecarInfo:
    annotation_id: str
    path: Path
    n_genes: int            # all genes in the shard
    n_genes_named: int      # genes written to #genes
    n_name_rows: int
    n_distinct_names: int
    n_casefold_changed: int
    n_skipped_names: int
    shard_fingerprint: str
    sidecar_fingerprint: str  # sha256 of the file bytes
    size_bytes: int
    meta: dict[str, str]


def sidecar_path_for(shard_path: Path) -> Path:
    return shard_path.parent / SIDECAR_NAME


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).translate(_BAD)


def _open_ro(shard_path: Path) -> sqlite3.Connection:
    uri = f"file:{shard_path}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.execute("PRAGMA query_only=ON")
    return conn


def shard_fingerprint(shard_path: Path, conn: sqlite3.Connection | None = None) -> str:
    """
    ``v1:<size>:<mtime_ns>:<n_genes>:<max_local_id>``.

    Cheap: one ``stat`` and two aggregate queries on the ``gene`` primary key.
    """
    st = shard_path.stat()
    own = conn is None
    if own:
        conn = _open_ro(shard_path)
    try:
        n, mx = conn.execute("SELECT COUNT(*), COALESCE(MAX(local_id), 0) FROM gene").fetchone()
    finally:
        if own:
            conn.close()
    return f"v1:{st.st_size}:{st.st_mtime_ns}:{n}:{mx}"


# Above this many name rows the dedupe/sort scratch table lives in a file (in
# ``tmp_dir``) instead of memory, so memory stays bounded for very large shards.
_MEMORY_ROWS = 1_000_000
_BATCH = 50_000


def _flush(scratch: sqlite3.Connection, batch: list[tuple[str, str, int]]) -> None:
    scratch.executemany("INSERT OR IGNORE INTO n VALUES (?, ?, ?)", batch)
    scratch.executemany("INSERT OR IGNORE INTO g VALUES (?)", [(b[2],) for b in batch])
    batch.clear()


def _build_sorted_names(
    conn: sqlite3.Connection, tmp_dir: Path
) -> tuple[sqlite3.Connection, Path | None, int, int, int, int]:
    """
    Fill a scratch SQLite table with distinct casefolded ``(name, kind, local_id)``.

    SQLite's BINARY collation compares UTF-8 bytes, which orders like Python's
    code-point order, so iterating the primary key yields the sidecar order without
    holding the rows in Python. Returns
    ``(scratch_conn, scratch_path, n_rows, n_distinct_names, n_changed, n_skipped)``.
    """
    n_src = conn.execute(
        "SELECT COUNT(*) FROM gene_xref WHERE namespace IN ('symbol', 'alias')"
    ).fetchone()[0]
    scratch_path: Path | None = None
    if n_src > _MEMORY_ROWS:
        tmp_dir.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(dir=tmp_dir, suffix=".names.sqlite")
        os.close(fd)
        scratch_path = Path(name)
        scratch = sqlite3.connect(name)
    else:
        scratch = sqlite3.connect(":memory:")
    scratch.execute("PRAGMA journal_mode=OFF")
    scratch.execute("PRAGMA synchronous=OFF")
    scratch.execute("PRAGMA cache_size=-65536")
    scratch.execute(
        "CREATE TABLE n (name TEXT NOT NULL, kind TEXT NOT NULL, local_id INTEGER NOT NULL, "
        "PRIMARY KEY (name, kind, local_id)) WITHOUT ROWID"
    )
    # Named gene ids, kept in a keyed table so later queries need no temporary sort or
    # DISTINCT structure (those can spill to the temp directory, possibly local disk).
    scratch.execute("CREATE TABLE g (local_id INTEGER PRIMARY KEY)")
    changed = 0
    skipped = 0
    batch: list[tuple[str, str, int]] = []
    cur = conn.execute(
        "SELECT accession, namespace, local_id FROM gene_xref "
        "WHERE namespace IN ('symbol', 'alias')"
    )
    for accession, namespace, local_id in cur:
        if accession is None:
            skipped += 1
            continue
        raw = str(accession)
        name = raw.casefold()
        if name != raw:
            changed += 1
        if not name or name != name.translate(_BAD):
            skipped += 1
            continue
        batch.append((name, _NAMESPACE_KIND[namespace], int(local_id)))
        if len(batch) >= _BATCH:
            _flush(scratch, batch)
    if batch:
        _flush(scratch, batch)
    scratch.commit()
    n_rows = scratch.execute("SELECT COUNT(*) FROM n").fetchone()[0]
    n_distinct = scratch.execute("SELECT COUNT(DISTINCT name) FROM n").fetchone()[0]
    return scratch, scratch_path, n_rows, n_distinct, changed, skipped


def _iter_named_genes(conn: sqlite3.Connection, scratch: sqlite3.Connection):
    """
    Yield gene rows, ordered by ``local_id``, only for ids that appear in the names
    scratch table. Uses primary-key lookups in bounded batches, so a shard with
    millions of unnamed genes is never scanned.
    """
    ids = scratch.execute("SELECT local_id FROM g ORDER BY local_id")
    while True:
        chunk = [r[0] for r in ids.fetchmany(500)]
        if not chunk:
            return
        marks = ",".join("?" for _ in chunk)
        yield from conn.execute(
            "SELECT local_id, source_gene_id, seqid, start, end, strand, feature_type, "
            f"biotype, primary_name FROM gene WHERE local_id IN ({marks}) ORDER BY local_id",
            chunk,
        )


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_sidecar(
    shard_path: Path,
    out_path: Path | None = None,
    *,
    built_at: str | None = None,
    tmp_dir: Path | None = None,
) -> SidecarInfo:
    """
    Write the sidecar for ``shard_path`` atomically and return its summary.

    Memory is bounded: names are deduplicated and sorted in a scratch SQLite table
    (in a file under ``tmp_dir`` for large shards, default: next to the output; keep
    it off local disk), genes and names are streamed into the gzip file.
    """
    shard_path = Path(shard_path)
    out_path = Path(out_path) if out_path else sidecar_path_for(shard_path)
    tmp_dir = Path(tmp_dir) if tmp_dir else out_path.parent
    built_at = built_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    conn = _open_ro(shard_path)
    scratch = None
    scratch_path: Path | None = None
    try:
        fingerprint = shard_fingerprint(shard_path, conn)
        meta = {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
        annotation_id = meta.get("annotation_id", "")
        if not annotation_id:
            raise SidecarError(f"shard has no annotation_id in meta: {shard_path}")
        n_genes = conn.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
        scratch, scratch_path, n_rows, n_distinct, changed, skipped = _build_sorted_names(
            conn, tmp_dir
        )
        n_named = scratch.execute("SELECT COUNT(*) FROM g").fetchone()[0]
        header = [
            ("sidecar_version", SIDECAR_VERSION),
            ("annotation_id", annotation_id),
            ("taxid", meta.get("taxid", "")),
            ("assembly_accession", meta.get("assembly_accession", "")),
            ("organism_name", _clean(meta.get("organism_name", ""))),
            ("source_database", _clean(meta.get("source_database", ""))),
            ("source_provider", _clean(meta.get("source_provider", ""))),
            ("gff_path", _clean(meta.get("gff_path", ""))),
            ("n_genes", n_genes),
        ("n_genes_named", n_named),
            ("n_name_rows", n_rows),
            ("n_distinct_names", n_distinct),
            ("n_casefold_changed", changed),
            ("n_skipped_names", skipped),
            ("shard_fingerprint", fingerprint),
            ("built_at", built_at),
        ]
        n_genes_written = 0
        with open(tmp, "wb") as fh:
            # mtime=0 and no filename keep the gzip bytes deterministic for equal content.
            with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=6) as gz:
                text = io.TextIOWrapper(gz, encoding="utf-8", newline="\n")
                for key, value in header:
                    text.write(f"#{key}={value}\n")
                text.write("#names\n")
                for name, kind, local_id in scratch.execute(
                    "SELECT name, kind, local_id FROM n ORDER BY name, kind, local_id"
                ):
                    text.write(f"{name}\t{kind}\t{local_id}\n")
                text.write("#genes\n")
                for g in _iter_named_genes(conn, scratch):
                    local_id, sgid, seqid, start, end, strand, ftype, biotype, pname = g
                    text.write(
                        f"{local_id}\t{_clean(sgid)}\t{_clean(seqid)}\t{start}\t{end}\t{strand}\t"
                        f"{_clean(ftype)}\t{_clean(biotype)}\t{_clean(pname)}\n"
                    )
                    n_genes_written += 1
                text.write(f"#end names={n_rows} genes={n_genes_written}\n")
                text.flush()
                text.detach()
            fh.flush()
            os.fsync(fh.fileno())
        if n_genes_written != n_named:
            raise SidecarError(
                f"shard changed while reading ({n_named} named genes expected, "
                f"{n_genes_written} written)"
            )
        os.replace(tmp, out_path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    finally:
        conn.close()
        if scratch is not None:
            scratch.close()
        if scratch_path is not None:
            scratch_path.unlink(missing_ok=True)

    return SidecarInfo(
        annotation_id=annotation_id,
        path=out_path,
        n_genes=n_genes,
        n_genes_named=n_genes_written,
        n_name_rows=n_rows,
        n_distinct_names=n_distinct,
        n_casefold_changed=changed,
        n_skipped_names=skipped,
        shard_fingerprint=fingerprint,
        sidecar_fingerprint=_sha256_file(out_path),
        size_bytes=out_path.stat().st_size,
        meta=meta,
    )


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def _lines(path: Path) -> Iterator[str]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="\n") as fh:
            for line in fh:
                yield line.rstrip("\n")
    except (OSError, EOFError) as exc:  # gzip errors, truncated stream
        raise SidecarError(f"cannot read sidecar {path}: {exc}") from exc


def read_header(path: Path) -> dict[str, str]:
    """Parse the ``#key=value`` header; stops at ``#names``."""
    header: dict[str, str] = {}
    for line in _lines(Path(path)):
        if line == "#names":
            break
        if not line.startswith("#") or "=" not in line:
            raise SidecarError(f"bad header line in {path}: {line[:60]!r}")
        key, _, value = line[1:].partition("=")
        header[key] = value
    else:
        raise SidecarError(f"no #names section in {path}")
    if header.get("sidecar_version") != str(SIDECAR_VERSION):
        raise SidecarError(
            f"unsupported sidecar_version {header.get('sidecar_version')!r} in {path}"
        )
    return header


def _iter_sections(path: Path) -> Iterator[tuple[str, list[str]]]:
    """Yield (section, fields) for every data row; verify the ``#end`` trailer."""
    section: str | None = None
    counts = {"names": 0, "genes": 0}
    ended = False
    in_header = True
    for line in _lines(Path(path)):
        if in_header:
            if line == "#names":
                in_header = False
                section = "names"
            continue
        if line == "#genes":
            section = "genes"
            continue
        if line.startswith("#end "):
            declared = dict(p.split("=", 1) for p in line[5:].split())
            if declared != {k: str(v) for k, v in counts.items()}:
                raise SidecarError(f"row counts differ from trailer in {path}")
            ended = True
            continue
        if ended:
            raise SidecarError(f"data after #end in {path}")
        if section is None:
            raise SidecarError(f"data before a section in {path}")
        counts[section] += 1
        yield section, line.split("\t")
    if not ended:
        raise SidecarError(f"sidecar truncated (no #end): {path}")


def iter_names(path: Path) -> Iterator[tuple[str, str, int]]:
    """Yield ``(name, kind, local_id)``; raises on a truncated or inconsistent file."""
    for section, f in _iter_sections(path):
        if section == "names":
            yield f[0], f[1], int(f[2])


def iter_genes(path: Path) -> Iterator[tuple]:
    """Yield gene rows as tuples (see module docstring); same integrity checks."""
    for section, f in _iter_sections(path):
        if section == "genes":
            yield (
                int(f[0]), f[1], f[2], int(f[3]), int(f[4]), int(f[5]),
                f[6], f[7], f[8],
            )


def is_stale(shard_path: Path, sidecar_path: Path | None = None) -> bool:
    """True when the sidecar is missing, unreadable or built from another shard state."""
    shard_path = Path(shard_path)
    sidecar_path = Path(sidecar_path) if sidecar_path else sidecar_path_for(shard_path)
    if not sidecar_path.is_file():
        return True
    try:
        header = read_header(sidecar_path)
    except SidecarError:
        return True
    return header.get("shard_fingerprint") != shard_fingerprint(shard_path)
