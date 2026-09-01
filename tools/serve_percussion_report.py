"""Serve a generated percussion report on every network interface."""

from __future__ import annotations

import argparse
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
import os


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--port", type=int, default=8000)
    arguments = parser.parse_args()
    directory = arguments.directory.resolve()
    if not directory.is_dir():
        parser.error(f"not a directory: {directory}")
    os.chdir(directory)
    server = ThreadingHTTPServer(("0.0.0.0", arguments.port), SimpleHTTPRequestHandler)
    print(f"serving {directory} at http://0.0.0.0:{arguments.port}/", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
