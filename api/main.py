"""Annotrieve-xrefs FastAPI application."""
from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from helpers.session import close_reader
from routes import router


@asynccontextmanager
async def lifespan(_app: FastAPI):
    yield
    close_reader()


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
        allow_methods=["GET", "OPTIONS", "HEAD"],
        allow_headers=["Content-Type", "Accept", "Cache-Control"],
        expose_headers=["Content-Length", "Content-Type", "Cache-Control"],
        max_age=86400,
    )

    app.include_router(router)
    return app


app = create_app()
