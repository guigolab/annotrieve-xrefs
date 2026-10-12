"""Map annotation.source_database to a compact integer for index payloads."""
from __future__ import annotations

# Matches gene_corpus.catalog.report.API_DB_SOURCES (case-insensitive).
SOURCE_NAME_TO_CODE: dict[str, int] = {
    "ensembl": 1,
    "genbank": 2,
    "refseq": 3,
}

# Legend keys are strings so JSON object keys stay uniform.
SOURCES_LEGEND: dict[str, str] = {
    "1": "Ensembl",
    "2": "GenBank",
    "3": "RefSeq",
}


def source_code(source_database: str | None) -> int:
    """Return 1/2/3 for Ensembl/GenBank/RefSeq, else 0."""
    if not source_database:
        return 0
    return SOURCE_NAME_TO_CODE.get(str(source_database).strip().casefold(), 0)
