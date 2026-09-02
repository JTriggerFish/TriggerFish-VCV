"""Sanitized browser index for a private local reference corpus."""

from __future__ import annotations

import json
from pathlib import Path

CORPUS_ID = "private-crash-a"


def load_corpus(root: Path) -> dict[str, object] | None:
    manifest = root / "cells.json"
    if not manifest.is_file():
        return None
    source = json.loads(manifest.read_text(encoding="utf-8"))
    cells = []
    for cell in source["cells"]:
        cells.append(
            {
                "label": cell["label"],
                "articulation": cell["articulation"],
                "velocity": cell["midi_velocity"],
                "repeat": cell["repeat"],
                "strength": cell["strength"],
                "location": cell["location"],
                "hardness": cell["hardness"],
                "seed": cell["seed"],
                "split": cell["split"],
                "url": f"/reference/{CORPUS_ID}/{cell['path']}",
            }
        )
    return {
        "id": CORPUS_ID,
        "name": "Private crash A",
        "audition_trim_db": 25.5,
        "cells": cells,
    }


def reference_path(root: Path, request_path: str) -> Path | None:
    prefix = f"/reference/{CORPUS_ID}/"
    if not request_path.startswith(prefix):
        return None
    name = request_path.removeprefix(prefix)
    if not name or Path(name).name != name:
        return None
    candidate = root / name
    return candidate if candidate.is_file() else None
