"""Turn deterministic cymbal captures into labeled crash fitting cells."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer, read_wav, write_wav
from .crash_fitting import CrashFitCell

ARTICULATION_CONTROLS = {
    "bell-tip": (0.0, 0.60),
    "bell-shank": (0.0, 0.92),
    "bow-tip": (0.5, 0.60),
    "bow-shank": (0.5, 0.92),
    "edge": (1.0, 0.60),
}


@dataclass(frozen=True)
class CapturedCell:
    label: str
    path: str
    articulation: str
    midi_velocity: int
    repeat: int
    strength: float
    location: float
    hardness: float
    seed: int
    split: str
    measured_onset_seconds: float
    cell_onset_seconds: float


def isolate_capture(
    audio: AudioBuffer,
    sweep_manifest: dict,
    output_directory: Path,
    duration_seconds: float = 10.0,
    search_seconds: float = 0.08,
    pre_roll_seconds: float = 0.05,
) -> tuple[CapturedCell, ...]:
    """Slice one lossless sweep render without normalizing individual hits."""
    mono = audio.mono()
    output_directory.mkdir(parents=True, exist_ok=True)
    sample_count = round(duration_seconds * mono.sample_rate)
    search = round(search_seconds * mono.sample_rate)
    pre_roll = round(pre_roll_seconds * mono.sample_rate)
    cells = []
    for hit in sweep_manifest["hits"]:
        nominal = round(float(hit["onset_seconds"]) * mono.sample_rate)
        search_start = max(0, nominal - search)
        search_end = min(mono.samples.size, nominal + search)
        local = detect_impact_onset(
            mono.samples[search_start:search_end], mono.sample_rate
        )
        onset = search_start + local
        segment_start = max(0, onset - pre_roll)
        destination_start = pre_roll - (onset - segment_start)
        segment = np.zeros(sample_count, dtype=np.float64)
        available = min(
            sample_count - destination_start, mono.samples.size - segment_start
        )
        if available > 0:
            segment[destination_start : destination_start + available] = mono.samples[
                segment_start : segment_start + available
            ]
        articulation = str(hit["articulation"])
        location, hardness = ARTICULATION_CONTROLS[articulation]
        velocity = int(hit["velocity"])
        repeat = int(hit["repeat"])
        name = (
            f"{int(hit['index']):03d}-{articulation}-v{velocity:03d}-r{repeat:02d}.wav"
        )
        write_wav(output_directory / name, AudioBuffer(segment, mono.sample_rate))
        cells.append(
            CapturedCell(
                label=f"{articulation} v{velocity:03d} r{repeat:02d}",
                path=name,
                articulation=articulation,
                midi_velocity=velocity,
                repeat=repeat,
                strength=velocity / 127.0,
                location=location,
                hardness=hardness,
                seed=0x53443300 + int(hit["index"]),
                split="fit" if repeat == 1 else "validation",
                measured_onset_seconds=onset / mono.sample_rate,
                cell_onset_seconds=pre_roll / mono.sample_rate,
            )
        )
    return tuple(cells)


def write_cells_manifest(
    path: Path, cells: tuple[CapturedCell, ...], source_audio: Path
) -> None:
    payload = {
        "schema": 1,
        "source_audio": str(source_audio.resolve()),
        "level_policy": "source gain preserved; cells are not normalized",
        "cells": [asdict(cell) for cell in cells],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def load_fit_cells(path: Path, split: str | None = None) -> tuple[CrashFitCell, ...]:
    """Load isolated cells while resolving WAV paths beside their manifest."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    cells = []
    for cell in payload["cells"]:
        if split is not None and cell["split"] != split:
            continue
        cells.append(
            CrashFitCell(
                cell["label"],
                read_wav(path.parent / cell["path"], "mean"),
                float(cell["strength"]),
                float(cell["location"]),
                float(cell["hardness"]),
                int(cell["seed"]),
            )
        )
    return tuple(cells)
