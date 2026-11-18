"""Genes a /sequences request refers to: from a ``query`` or from explicit ``locus_ref`` entries.

Both paths read only ``presence.sqlite``, ``loci.sqlite`` and ``loci_ref.sqlite``. No old shard or
corpus table is touched.
"""
from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass

from helpers.presence import annotation_info, annotation_paths
from services.presence_loci import STRAND, LociQuery, build_query, estimate, iter_loci
from services.seq_request import LociEntry, SeqRequest


@dataclass
class Gene:
    ord: int
    annotation_id: str
    local_id: int
    organism: str
    accession: str
    taxid: int
    gff_path: str
    seqid: str
    start: int
    end: int
    strand: str
    symbol: str
    gene_id: str
    pinned: str | None = None

    @property
    def locus_ref(self) -> str:
        return f"{self.annotation_id}:{self.local_id}"


def build_loci_query(pres: sqlite3.Connection, reader, data_version: str, req: SeqRequest
                     ) -> LociQuery:
    q = req.query or {}
    return build_query(
        pres, reader.all_mask, data_version, q["names"], terms=q.get("terms") or None,
        include_taxids=q["include_taxids"], exclude_taxids=q["exclude_taxids"], alias=q["alias"],
        annotation_ids=q["annotation_ids"], case_sensitive=q["case_sensitive"], flank_bp=0)


def query_estimate(pres, loci, q: LociQuery) -> dict:
    return estimate(pres, loci, q)


def query_genes(pres, loci, reader, q: LociQuery) -> Iterator[Gene]:
    info = annotation_info(reader)
    paths = annotation_paths(reader)
    for row in iter_loci(pres, loci, q, None, info):
        taxid, gff_path = paths[row["ord"]]
        yield Gene(
            ord=row["ord"], annotation_id=row["annotation_id"], local_id=row["local_id"],
            organism=row["organism_name"] or "", accession=row["assembly_accession"] or "",
            taxid=taxid, gff_path=gff_path, seqid=row["seqid"], start=row["start"],
            end=row["end"], strand=row["strand"], symbol=row["symbol"] or "",
            gene_id=row["gene_id"] or "")


def ref_genes(pres, loci, ref, reader, entries: tuple[LociEntry, ...]
              ) -> Iterator[tuple[LociEntry, Gene | None, str]]:
    """For each entry: ``(entry, gene, "")`` or ``(entry, None, reason)``.

    Reasons: ``annotation_not_found`` (not in this index) and ``locus_not_found``.
    """
    info = annotation_info(reader)
    paths = annotation_paths(reader)
    ords: dict[str, int | None] = {}
    seqids: dict[tuple[int, int], str] = {}
    for entry in entries:
        if entry.annotation_id not in ords:
            row = pres.execute("SELECT ord FROM annotation WHERE annotation_id=?",
                               (entry.annotation_id,)).fetchone()
            ords[entry.annotation_id] = int(row[0]) if row else None
        ord_ = ords[entry.annotation_id]
        if ord_ is None:
            yield entry, None, "annotation_not_found"
            continue
        ref_row = ref.execute("SELECT name_id FROM ref WHERE ord=? AND local_id=?",
                              (ord_, entry.local_id)).fetchone()
        r = None if ref_row is None else loci.execute(
            "SELECT seqid_id, start, end, strand, gene_id, primary_name FROM locus "
            "WHERE name_id=? AND ord=? AND local_id=?",
            (ref_row[0], ord_, entry.local_id)).fetchone()
        if r is None:
            yield entry, None, "locus_not_found"
            continue
        sq = (ord_, r["seqid_id"])
        if sq not in seqids:
            row = loci.execute("SELECT seqid FROM seqid WHERE ord=? AND seqid_id=?", sq).fetchone()
            seqids[sq] = row[0] if row else ""
        _aid, organism, accession, _source = info[ord_]
        taxid, gff_path = paths[ord_]
        gene_id = r["gene_id"] or ""
        yield entry, Gene(
            ord=ord_, annotation_id=entry.annotation_id, local_id=entry.local_id,
            organism=organism or "", accession=accession or "", taxid=taxid, gff_path=gff_path,
            seqid=seqids[sq], start=r["start"], end=r["end"],
            strand=STRAND.get(r["strand"], "."), symbol=r["primary_name"] or gene_id,
            gene_id=gene_id, pinned=entry.transcript_id), ""
