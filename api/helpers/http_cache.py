"""ETag and Cache-Control for answers that depend only on the published data version."""
from __future__ import annotations

from starlette.concurrency import run_in_threadpool

from helpers.presence import current_data_version, data_version_check_due, refresh_data_version
from settings import settings

# Routes whose answer is a pure function of the URL and the data version.
_PREFIXES = ("/hits", "/names", "/loci", "/distributions", "/annotations", "/genes", "/coverage")


class DataVersionCache:
    """
    Pure ASGI middleware. For GET/HEAD on the data routes:

    * a request whose ``If-None-Match`` carries the current version gets ``304`` at once, with no
      query run;
    * a 200 answer gets ``ETag: W/"<version>"`` and ``Cache-Control: public, max-age=N`` (N from
      ``HTTP_CACHE_MAX_AGE_S``), unless the data version changed while it was being produced, in
      which case it gets ``no-store``.
    """

    def __init__(self, app) -> None:
        self.app = app

    async def __call__(self, scope, receive, send):
        if (scope["type"] != "http" or scope["method"] not in ("GET", "HEAD")
                or not scope["path"].startswith(_PREFIXES)):
            return await self.app(scope, receive, send)
        if data_version_check_due():
            # File access (network drive): off the event loop, and only every few seconds.
            await run_in_threadpool(refresh_data_version)
        version = current_data_version()
        if not version:
            return await self.app(scope, receive, send)
        etag = f'W/"{version}"'
        cache = (f"public, max-age={settings.HTTP_CACHE_MAX_AGE_S}"
                 if settings.HTTP_CACHE_MAX_AGE_S > 0 else None)
        sent = dict(scope.get("headers") or [])
        if etag.encode() in [p.strip() for p in sent.get(b"if-none-match", b"").split(b",")]:
            headers = [(b"etag", etag.encode())]
            if cache:
                headers.append((b"cache-control", cache.encode()))
            await send({"type": "http.response.start", "status": 304, "headers": headers})
            await send({"type": "http.response.body", "body": b""})
            return

        async def send_with_headers(message):
            if message["type"] == "http.response.start" and message["status"] == 200:
                headers = list(message.get("headers") or [])
                names = {k.lower() for k, _ in headers}
                if current_data_version() != version:
                    headers.append((b"cache-control", b"no-store"))
                else:
                    if b"etag" not in names:
                        headers.append((b"etag", etag.encode()))
                    if cache and b"cache-control" not in names:
                        headers.append((b"cache-control", cache.encode()))
                message = {**message, "headers": headers}
            await send(message)

        await self.app(scope, receive, send_with_headers)
