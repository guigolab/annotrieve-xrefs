"""Write the API's OpenAPI schema to ``api/openapi.json`` (or print it with --stdout).

    python3 tools/export_openapi.py

A test fails when the committed file differs from the app, so run this after changing a route.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "api"), str(ROOT)]


def render() -> str:
    from main import app
    return json.dumps(app.openapi(), indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    if "--stdout" in sys.argv:
        sys.stdout.write(render())
    else:
        (ROOT / "api" / "openapi.json").write_text(render())
        print("wrote api/openapi.json")
