"""Serve the optional percussion workbench with worker-isolation headers."""

from __future__ import annotations

import argparse
import json
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Lock
from urllib.parse import urlsplit

from percussion_reference_corpus import build_catalog, reference_path


class ReferenceCatalog:
    """Atomically refresh the development-only reference allow-list."""

    def __init__(
        self,
        private_root: Path,
        library_root: Path | None,
        reference_root: Path | None,
    ) -> None:
        self._roots = (private_root, library_root, reference_root)
        self._lock = Lock()
        self._corpora: list[dict[str, object]] = []
        self._paths: dict[str, Path] = {}
        self.refresh()

    def refresh(self) -> list[dict[str, object]]:
        corpora, paths = build_catalog(*self._roots)
        with self._lock:
            self._corpora, self._paths = corpora, paths
            return self._corpora

    def resolve(self, request_path: str) -> Path | None:
        with self._lock:
            return reference_path(self._paths, request_path)


class WorkbenchHandler(SimpleHTTPRequestHandler):
    def __init__(
        self,
        *args,
        catalog: ReferenceCatalog,
        **kwargs,
    ) -> None:
        self.catalog = catalog
        super().__init__(*args, **kwargs)

    def do_GET(self) -> None:
        if urlsplit(self.path).path == "/api/reference-corpora":
            self._send_json({"corpora": self.catalog.refresh()})
            return
        super().do_GET()

    def translate_path(self, path: str) -> str:
        reference = self.catalog.resolve(urlsplit(path).path)
        return str(reference) if reference else super().translate_path(path)

    def _send_json(self, value: object) -> None:
        payload = json.dumps(value, separators=(",", ":")).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def end_headers(self) -> None:
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument(
        "--corpus-root",
        type=Path,
        default=Path(
            "build/cymbal-calibration/references/"
            "private-corpus-a-crash-v1/cells-oh-dyn-v2"
        ),
    )
    parser.add_argument("--library-root", type=Path)
    parser.add_argument("--reference-root", type=Path)
    arguments = parser.parse_args()
    directory = arguments.directory.resolve()
    if not directory.is_dir():
        parser.error(f"not a directory: {directory}")
    corpus_root = arguments.corpus_root.resolve()
    catalog = ReferenceCatalog(
        corpus_root, arguments.library_root, arguments.reference_root
    )
    handler = partial(
        WorkbenchHandler,
        directory=str(directory),
        catalog=catalog,
    )
    server = ThreadingHTTPServer(("0.0.0.0", arguments.port), handler)
    print(f"serving {directory} at http://0.0.0.0:{arguments.port}/", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
