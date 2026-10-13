"""Aggregate API routers (nginx strips /annotrieve-xrefs/api/v1).

Do not use ``from __future__ import annotations`` here: it binds the name
``annotations`` and breaks ``from routers import annotations``.
"""

from fastapi import APIRouter

from routers import annotations, health, namespaces

api_router = APIRouter()
api_router.include_router(health.router)
# /matches before namespaced paths that could shadow (registered on namespaces).
api_router.include_router(namespaces.router)
api_router.include_router(annotations.router)
