"""Aggregate API routers (nginx strips /annotrieve-xrefs/api/v1).

Do not use ``from __future__ import annotations`` here: it binds the name
``annotations`` and breaks ``from routers import annotations``.
"""

from fastapi import APIRouter

from routers import annotations, coverage, distributions, downloads, health, hits, loci, names

api_router = APIRouter()
api_router.include_router(health.router, tags=["system"])
api_router.include_router(hits.router, tags=["search"])
api_router.include_router(names.router, tags=["search"])
api_router.include_router(coverage.router, tags=["search"])
api_router.include_router(loci.router, tags=["loci"])
api_router.include_router(distributions.router, tags=["distributions"])
api_router.include_router(annotations.router, tags=["annotations"])
api_router.include_router(downloads.router, tags=["downloads"])
