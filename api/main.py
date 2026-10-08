"""Annotrieve-xrefs FastAPI application."""
from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from helpers.errors import error_detail
from helpers.session import close_reader
from routers import api_router


@asynccontextmanager
async def lifespan(_app: FastAPI):
    yield
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
    app = FastAPI(
        title="Annotrieve xrefs",
        description=(
            "Search Annotrieve annotation files by gene symbol, gene id, "
            "or other accession, and see the matching genes."
        ),
        version="1.0.0",
        lifespan=lifespan,
    )

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
        expose_headers=["Content-Length", "Content-Type", "Cache-Control"],
        max_age=86400,
    )

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
