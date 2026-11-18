"""The client's IP behind the proxy chain (Traefik -> nginx -> API)."""
from __future__ import annotations

import hashlib

from fastapi import Request

from settings import settings


def client_ip(request: Request) -> str:
    """Address of the caller.

    ``X-Forwarded-For`` is a list the proxies append to; a client may send its own leading entries,
    so we take the one ``TRUSTED_PROXY_HOPS`` places from the right (the address the first trusted
    proxy saw). With no header, or hops=0, the socket peer is used.
    """
    hops = int(settings.TRUSTED_PROXY_HOPS)
    if hops > 0:
        parts = [p.strip() for p in request.headers.get("x-forwarded-for", "").split(",") if p.strip()]
        if len(parts) >= hops:
            return parts[-hops]
    return request.client.host if request.client else "unknown"


def owner_key(ip: str) -> str:
    """What is stored next to a result: a hash, not the address."""
    return hashlib.sha256(("annotrieve-xrefs|" + ip).encode()).hexdigest()[:16]
