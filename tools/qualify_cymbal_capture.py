"""Qualify an isolated cymbal grid before admitting it to fitting."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from triggerfish_percussion.capture_quality import qualify_cells


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = qualify_cells(arguments.cells_manifest)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))
    print("ACCEPTED" if result["accepted"] else "REJECTED")
    if not result["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
