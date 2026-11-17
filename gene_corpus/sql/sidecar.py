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
* ``#genes`` is sorted by ``local_id``. NULL text values are written as an empty field.
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
    n_genes: int
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


def _collect_names(conn: sqlite3.Connection) -> tuple[list[tuple[str, str, int]], int, int]:
    """Return sorted distinct (name, kind, local_id), n_casefold_changed, n_skipped."""
    seen: set[tuple[str, str, int]] = set()
    changed = 0
    skipped = 0
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
        seen.add((name, _NAMESPACE_KIND[namespace], int(local_id)))
    return sorted(seen), changed, skipped


def write_sidecar(
    shard_path: Path,
    out_path: Path | None = None,
    *,
    built_at: str | None = None,
) -> SidecarInfo:
    """Write the sidecar for ``shard_path`` atomically and return its summary."""
    shard_path = Path(shard_path)
    out_path = Path(out_path) if out_path else sidecar_path_for(shard_path)
    built_at = built_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    conn = _open_ro(shard_path)
    try:
        fingerprint = shard_fingerprint(shard_path, conn)
        meta = {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
        annotation_id = meta.get("annotation_id", "")
        if not annotation_id:
            raise SidecarError(f"shard has no annotation_id in meta: {shard_path}")
        names, changed, skipped = _collect_names(conn)
        genes = conn.execute(
            "SELECT local_id, source_gene_id, seqid, start, end, strand, feature_type, "
            "biotype, primary_name FROM gene ORDER BY local_id"
        ).fetchall()
    finally:
        conn.close()

    n_distinct = len({n for n, _k, _l in names})
    header = [
        ("sidecar_version", SIDECAR_VERSION),
        ("annotation_id", annotation_id),
        ("taxid", meta.get("taxid", "")),
        ("assembly_accession", meta.get("assembly_accession", "")),
        ("organism_name", _clean(meta.get("organism_name", ""))),
        ("source_database", _clean(meta.get("source_database", ""))),
        ("source_provider", _clean(meta.get("source_provider", ""))),
        ("gff_path", _clean(meta.get("gff_path", ""))),
        ("n_genes", len(genes)),
        ("n_name_rows", len(names)),
        ("n_distinct_names", n_distinct),
        ("n_casefold_changed", changed),
        ("n_skipped_names", skipped),
        ("shard_fingerprint", fingerprint),
        ("built_at", built_at),
    ]

    buf = io.BytesIO()
    # mtime=0 and no filename keep the gzip bytes deterministic for equal content.
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0, compresslevel=6) as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="\n")
        for key, value in header:
            text.write(f"#{key}={value}\n")
        text.write("#names\n")
        for name, kind, local_id in names:
            text.write(f"{name}\t{kind}\t{local_id}\n")
        text.write("#genes\n")
        for g in genes:
            local_id, sgid, seqid, start, end, strand, ftype, biotype, pname = g
            text.write(
                f"{local_id}\t{_clean(sgid)}\t{_clean(seqid)}\t{start}\t{end}\t{strand}\t"
                f"{_clean(ftype)}\t{_clean(biotype)}\t{_clean(pname)}\n"
            )
        text.write(f"#end names={len(names)} genes={len(genes)}\n")
        text.flush()
        text.detach()
    data = buf.getvalue()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, out_path)

    return SidecarInfo(
        annotation_id=annotation_id,
        path=out_path,
        n_genes=len(genes),
        n_name_rows=len(names),
        n_distinct_names=n_distinct,
        n_casefold_changed=changed,
        n_skipped_names=skipped,
        shard_fingerprint=fingerprint,
        sidecar_fingerprint=hashlib.sha256(data).hexdigest(),
        size_bytes=len(data),
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
