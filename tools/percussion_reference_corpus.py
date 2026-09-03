"""Build a small, sanitized index over local percussion references."""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.parse import quote, unquote

PRIVATE_CRASH_ID = "private-crash-a"


def _cell(
    corpus_id: str,
    path: Path,
    label: str,
    articulation: str,
    velocity: int,
    strength: float,
    location: float = 0.5,
    repeat: int = 1,
    hardness: float = 0.65,
    implement: float = 1.0,
    contact_spread: float = 0.2,
) -> dict[str, object]:
    return {
        "label": label,
        "articulation": articulation,
        "velocity": velocity,
        "repeat": repeat,
        "strength": strength,
        "location": location,
        "hardness": hardness,
        "implement": implement,
        "contactSpread": contact_spread,
        "onset_seconds": 0.0,
        "seed": 1000 + velocity * 7 + repeat,
        "split": "audition",
        "url": f"/reference/{corpus_id}/{quote(path.name)}",
    }


def _private_crash(root: Path) -> tuple[dict[str, object], dict[str, Path]] | None:
    manifest = root / "cells.json"
    if not manifest.is_file():
        return None
    source = json.loads(manifest.read_text(encoding="utf-8"))
    cells = []
    paths = {}
    selected_velocities = {24, 48, 72, 96, 120}
    for source_cell in source["cells"]:
        if (
            source_cell["midi_velocity"] not in selected_velocities
            or source_cell["repeat"] != 1
        ):
            continue
        path = (root / source_cell["path"]).resolve()
        if not path.is_file() or root.resolve() not in path.parents:
            continue
        cell = _cell(
            PRIVATE_CRASH_ID,
            path,
            source_cell["label"],
            source_cell["articulation"],
            source_cell["midi_velocity"],
            source_cell["strength"],
            source_cell["location"],
            source_cell["repeat"],
            source_cell["hardness"],
        )
        cell["onset_seconds"] = source_cell["cell_onset_seconds"]
        cell["seed"] = source_cell["seed"]
        cell["split"] = source_cell["split"]
        cells.append(cell)
        paths[unquote(cell["url"])] = path
    return (
        {
            "id": PRIVATE_CRASH_ID,
            "name": "Private crash A",
            "audition_trim_db": 25.5,
            "calibration": {
                "id": "crash-standard",
                "name": "Crash — medium edge",
                "recipe": "metal.cymbal.v1",
                "parameter_preset": "metallic-default",
                "articulation": "edge",
                "velocity": 72,
                "repeat": 1,
            },
            "cells": cells,
        },
        paths,
    )


def _curated(
    corpus_id: str,
    name: str,
    specifications: list[tuple[Path, str, str, int, float, float, int]],
    trim_db: float = 0.0,
    calibration: dict[str, object] | None = None,
    hardness: float = 0.65,
    implement: float = 1.0,
    contact_spread: float = 0.2,
) -> tuple[dict[str, object], dict[str, Path]] | None:
    cells = []
    paths = {}
    for (
        path,
        label,
        articulation,
        velocity,
        strength,
        location,
        repeat,
    ) in specifications:
        if not path.is_file():
            continue
        cell = _cell(
            corpus_id,
            path,
            label,
            articulation,
            velocity,
            strength,
            location,
            repeat,
            hardness,
            implement,
            contact_spread,
        )
        cells.append(cell)
        paths[unquote(cell["url"])] = path.resolve()
    if not cells:
        return None
    return (
        {
            "id": corpus_id,
            "name": name,
            "audition_trim_db": trim_db,
            "calibration": calibration,
            "cells": cells,
        },
        paths,
    )


def _installed_library_root() -> Path:
    return Path(os.environ.get("LOCALAPPDATA", "")) / (
        "Bitwig Studio/installed-packages/5.0/Bitwig/" "Acoustic Drums and Percussion"
    )


def _snare(library: Path):
    root = library / "Acoustic Drums/Snares/Yamaha Mapple CA"
    layers = [("A", 40, 0.32), ("B", 82, 0.65), ("C", 116, 0.92)]
    # The package identifies A/B/C as distinct dynamic layers. Its numbered
    # groups are not documented as strike positions, so do not invent a
    # location mapping for them: expose one consistent group until metadata is
    # available.
    positions = [(1, "main", 0.3)]
    specs = []
    for position, articulation, location in positions:
        for layer, velocity, strength in layers:
            filename = f"Snare YMCA {position} {layer} 01.wav"
            specs.append(
                (
                    root / filename,
                    filename.removesuffix(".wav"),
                    articulation,
                    velocity,
                    strength,
                    location,
                    1,
                )
            )
    return _curated(
        "acoustic-snare-maple",
        "Acoustic snare — maple",
        specs,
        5.0,
        {
            "id": "snare-standard",
            "name": "Snare — medium centre",
            "recipe": "drum.snare.v1",
            "parameter_preset": "snare-default",
            "articulation": "main",
            "velocity": 82,
            "repeat": 1,
        },
    )


def _kick(library: Path):
    root = library / "Acoustic Drums/Kicks/Yamaha Oak Custom"
    layers = [(1, 32, 0.25), (4, 64, 0.5), (7, 96, 0.76), (10, 122, 0.97)]
    specs = []
    for layer, velocity, strength in layers:
        filename = f"Kick Yamaha Oak Custom {layer} 01.wav"
        specs.append(
            (
                root / filename,
                filename.removesuffix(".wav"),
                "centre",
                velocity,
                strength,
                0.5,
                1,
            )
        )
    return _curated(
        "acoustic-kick-oak",
        "Acoustic kick — oak",
        specs,
        2.0,
        {
            "id": "kick-standard",
            "name": "Kick — medium centre",
            "recipe": "drum.kick-fm.v1",
            "parameter_preset": "kick-default",
            "articulation": "centre",
            "velocity": 64,
            "repeat": 1,
        },
        hardness=0.5,
        implement=0.5,
    )


def _gong(library: Path):
    root = library / "Percussion/Gongs"
    specs = []
    for repeat, number in enumerate((1, 2, 3, 5, 7), start=1):
        filename = f"Gong Dresden {number:02d}.wav"
        specs.append(
            (
                root / filename,
                filename.removesuffix(".wav"),
                "mallet",
                96,
                0.76,
                0.55,
                repeat,
            )
        )
    return _curated(
        "gong-dresden",
        "Dresden gong",
        specs,
        -6.0,
        {
            "id": "gong-standard",
            "name": "Gong — representative mallet",
            "recipe": "metal.cymbal.v1",
            "parameter_preset": "gong-start",
            "articulation": "mallet",
            "velocity": 96,
            "repeat": 3,
        },
        hardness=0.35,
        implement=0.5,
        contact_spread=0.3,
    )


def _ride(reference_root: Path):
    layers = [("pp", 40, 0.32), ("mf", 82, 0.65), ("ff", 116, 0.92)]
    positions = [
        ("bell", "bell", 0.08),
        ("normal", "bow", 0.58),
        ("shoulder", "shoulder", 0.86),
    ]
    specs = []
    for token, articulation, location in positions:
        for dynamic, velocity, strength in layers:
            filename = f"21ride.stick.{token}.{dynamic}.stereo.wav"
            specs.append(
                (
                    reference_root / filename,
                    filename.removesuffix(".wav"),
                    articulation,
                    velocity,
                    strength,
                    location,
                    1,
                )
            )
    return _curated(
        "ride-21-reference",
        "21-inch ride reference",
        specs,
        -3.0,
        {
            "id": "ride-standard",
            "name": "Ride — medium bow",
            "recipe": "metal.cymbal.v1",
            "parameter_preset": "ride-start",
            "articulation": "bow",
            "velocity": 82,
            "repeat": 1,
        },
    )


def build_catalog(
    private_root: Path,
    library_root: Path | None = None,
    reference_root: Path | None = None,
) -> tuple[list[dict[str, object]], dict[str, Path]]:
    """Return public corpus metadata and an exact allow-list of audio paths."""
    library = (library_root or _installed_library_root()).resolve()
    references = (reference_root or private_root.parents[1]).resolve()
    candidates = [
        _private_crash(private_root.resolve()),
        _snare(library),
        _kick(library),
        _gong(library),
        _ride(references),
    ]
    corpora = []
    paths: dict[str, Path] = {}
    for candidate in candidates:
        if candidate is None:
            continue
        corpus, corpus_paths = candidate
        corpora.append(corpus)
        paths.update(corpus_paths)
    return corpora, paths


def reference_path(paths: dict[str, Path], request_path: str) -> Path | None:
    """Resolve only paths emitted by :func:`build_catalog`."""
    return paths.get(unquote(request_path))
