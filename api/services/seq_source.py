"""Open the index connections for a /sequences request and produce its genes."""
from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass

from helpers.cursor import QueryError
from helpers.presence import (
    open_ref_connection,
    open_request_connections,
    presence_enabled,
    presence_reader,
)
from services.seq_export import Item
from services.seq_genes import build_loci_query, query_genes, ref_genes
from services.seq_request import SeqRequest
from services.seq_resolve import SeqError


def _convert(exc: Exception) -> SeqError:
    """The helpers raise FastAPI HTTPExceptions; services stay framework-free."""
    detail = getattr(exc, "detail", None)
    code = detail.get("code", "index_unavailable") if isinstance(detail, dict) else "index_unavailable"
    message = detail.get("message", str(exc)) if isinstance(detail, dict) else str(exc)
    return SeqError(code, message, status=getattr(exc, "status_code", 503))


@dataclass
class Counted:
    n_genes: int
    complete: bool           # False when counting stopped at the cap
    bp_upper_bound: int
    n_annotations: int | None = None


class Source:
    """Connections + gene stream of one request. Always call :meth:`close`."""

    def __init__(self, req: SeqRequest) -> None:
        if not presence_enabled():
            raise SeqError("sequences_unavailable", "The presence/loci index is not enabled",
                           status=503)
        self.req = req
        self.pres = self.loci = self.ref = None
        try:
            self.reader = presence_reader()
            self.pres, self.loci, self.data_version = open_request_connections()
            self.q = None
            if req.query is not None:
                try:
                    self.q = build_loci_query(self.pres, self.reader, self.data_version, req)
                except QueryError:
                    raise
            else:
                self.ref = open_ref_connection()
        except SeqError:
            self.close()
            raise
        except QueryError:
            self.close()
            raise
        except Exception as exc:                          # HTTPException from the helpers
            self.close()
            raise _convert(exc) from exc

    def close(self) -> None:
        for conn in (self.pres, self.loci, self.ref):
            if isinstance(conn, sqlite3.Connection):
                conn.close()
        self.pres = self.loci = self.ref = None

    def terms(self) -> list[dict]:
        if self.q is None:
            return []
        return [{"index": i, "input": t["input"], "name": t["name"],
                 "found": self.q.name_ids[i] is not None} for i, t in enumerate(self.q.terms)]

    def items(self) -> Iterator[Item]:
        if self.q is not None:
            for gene in query_genes(self.pres, self.loci, self.reader, self.q):
                yield gene.locus_ref, gene, ""
        else:
            for entry, gene, reason in ref_genes(self.pres, self.loci, self.ref, self.reader,
                                                 self.req.loci):
                yield entry.locus_ref, gene, reason

    def count(self, cap: int) -> Counted:
        """Count genes up to ``cap`` (+1 to detect overflow); cheap compared with the export."""
        n = bp = 0
        for _ref, gene, _reason in self.items():
            n += 1
            if gene is not None:
                bp += max(0, gene.end - gene.start + 1)
            if n > cap:
                return Counted(n, False, bp)
        return Counted(n, True, bp)
