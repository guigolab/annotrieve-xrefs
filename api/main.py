"""Annotrieve-xrefs FastAPI application."""
from __future__ import annotations

import os
import sqlite3
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from helpers.errors import error_detail
from helpers.http_cache import DataVersionCache
from helpers.logging_config import configure_logging
from helpers.ncbi_eutils import close_ncbi_client
from helpers.session import close_reader
from routers import api_router
from settings import settings

# Before create_app / worker serve: INFO must reach stderr for ``docker logs``.
configure_logging()


@asynccontextmanager
async def lifespan(_app: FastAPI):
    # Re-run after worker boot — gunicorn may have reconfigured logging.
    configure_logging()
    import anyio.to_thread
    anyio.to_thread.current_default_thread_limiter().total_tokens = max(4, settings.THREADPOOL_SIZE)
    # Scratch directory on the network drive (compose sets TMPDIR/SQLITE_TMPDIR there): make sure it exists.
    for var in ("TMPDIR", "SQLITE_TMPDIR"):
        if os.getenv(var):
            try:
                os.makedirs(os.environ[var], exist_ok=True)
            except OSError:
                pass
    # Open the published index now, so its data version is known before the first request.
    from starlette.concurrency import run_in_threadpool
    from helpers.presence import presence_enabled, presence_reader
    if presence_enabled():
        try:
            await run_in_threadpool(presence_reader)
        except Exception:  # not published yet: /health reports it
            pass
    yield
    close_ncbi_client()
    close_reader()


def _validation_message(errors: list[dict]) -> str:
    """Short human-readable summary from the first validation error."""
    if not errors:
        return "Invalid request"
    first = errors[0]
    loc = ".".join(str(part) for part in first.get("loc", ()))
    msg = first.get("msg") or "Invalid request"
    if loc:
        return f"{msg}: {loc}"
    return msg


def _compact_validation_errors(exc: RequestValidationError) -> list[dict]:
    """Keep loc / type / msg only (JSON-safe, no pydantic ctx objects)."""
    compact: list[dict] = []
    for err in exc.errors():
        compact.append(
            {
                "loc": list(err.get("loc", ())),
                "type": err.get("type"),
                "msg": err.get("msg"),
            }
        )
    return compact


def create_app() -> FastAPI:
    configure_logging()
    app = FastAPI(
        title="Annotrieve xrefs",
        description=(
            "From gene names to annotations, gene coordinates and sequences.\n\n"
            "- **search**: `/names` (prefix discovery), `/hits` (which annotations contain "
            "these exact gene names).\n"
            "- **loci**: `/loci` (where those genes are: one row per gene, as json, tsv "
            "or ndjson).\n"
            "- **sequences**: `/sequences` (CDS, transcript or protein FASTA for a query or a "
            "list of `locus_ref`s), with `/sequence-estimates`, `/sequence-plans`, jobs and "
            "cached results.\n"
            "- **annotations**: per-annotation genes and cross-references (legacy shard-backed)."
        ),
        openapi_tags=[
            {"name": "system", "description": "Liveness."},
            {"name": "search", "description": "Gene-name discovery and annotation presence."},
            {"name": "loci", "description": "Gene coordinates."},
            {"name": "sequences", "description": "FASTA export, previews, jobs and results."},
            {"name": "annotations", "description": "Per-annotation genes and xrefs (legacy)."},
        ],
        version="1.1.0",
        lifespan=lifespan,
    )

    app.add_middleware(DataVersionCache)

    allowed_origins_env = os.getenv("CORS_ALLOWED_ORIGINS", "")
    if allowed_origins_env:
        allowed_origins = [
            origin.strip()
            for origin in allowed_origins_env.split(",")
            if origin.strip()
        ]
    else:
        allowed_origins = ["*"]

    app.add_middleware(
        CORSMiddleware,
        allow_origins=allowed_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST", "OPTIONS", "HEAD"],
        allow_headers=["Content-Type", "Accept", "Cache-Control"],
        expose_headers=[
            "Content-Length",
            "Content-Type",
            "Cache-Control",
            "Content-Disposition",
            "X-Sequence-Count",
            "X-Skipped-Count",
            "X-Gene-Count",
            "X-Result-Id",
            "X-Data-Version",
            "X-Manifest-Url",
            "Location",
            "Retry-After",
        ],
        max_age=86400,
    )

    @app.exception_handler(sqlite3.OperationalError)
    async def sqlite_error_handler(_request: Request, exc: sqlite3.OperationalError) -> JSONResponse:
        if "interrupted" in str(exc):  # query time budget (services.annotations.open_shard)
            return JSONResponse(
                status_code=503,
                content={"detail": error_detail(
                    "query_too_slow", "the query took too long; narrow it or retry later")},
                headers={"Retry-After": "5"},
            )
        return JSONResponse(status_code=500, content={"detail": "Internal Server Error"})

    @app.exception_handler(RequestValidationError)
    async def validation_error_handler(
        _request: Request, exc: RequestValidationError
    ) -> JSONResponse:
        errors = _compact_validation_errors(exc)
        return JSONResponse(
            status_code=422,
            content={
                "detail": error_detail(
                    "validation_error",
                    _validation_message(errors),
                    errors=errors,
                )
            },
        )

    app.include_router(api_router)
    return app


app = create_app()
