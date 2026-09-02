"""Qualification measurements for locally rendered cymbal capture grids."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np

from .alignment import detect_impact_onset
from .audio_io import read_wav


@dataclass(frozen=True)
class CellQuality:
    articulation: str
    velocity: int
    repeat: int
    peak: float
    rms_100ms: float
    centroid_hz: float
    rolloff_85_hz: float
    onset_error_ms: float
    content_hash: str
    band_profile_db: tuple[float, ...]


def measure_cell(path: Path, metadata: dict) -> CellQuality:
    """Measure one level-preserved, onset-aligned capture cell."""
    audio = read_wav(path, "mean")
    expected = round(float(metadata["cell_onset_seconds"]) * audio.sample_rate)
    measured = detect_impact_onset(audio.samples, audio.sample_rate)
    count = round(0.100 * audio.sample_rate)
    signal = _fixed_window(audio.samples, expected, count)
    power = np.square(np.abs(np.fft.rfft(signal * np.hanning(signal.size))))
    frequencies = np.fft.rfftfreq(signal.size, 1.0 / audio.sample_rate)
    total = max(float(np.sum(power)), np.finfo(float).tiny)
    cumulative = np.cumsum(power)
    rolloff = frequencies[
        min(np.searchsorted(cumulative, 0.85 * total), frequencies.size - 1)
    ]
    profile = _band_profile(frequencies, power)
    digest = sha256(np.ascontiguousarray(audio.samples).view(np.uint8)).hexdigest()
    return CellQuality(
        str(metadata["articulation"]),
        int(metadata["midi_velocity"]),
        int(metadata["repeat"]),
        float(np.max(np.abs(signal))),
        float(np.sqrt(np.mean(np.square(signal)))),
        float(np.sum(frequencies * power) / total),
        float(rolloff),
        1000.0 * (measured - expected) / audio.sample_rate,
        digest,
        tuple(float(value) for value in profile),
    )


def qualify_cells(manifest_path: Path) -> dict:
    """Return reproducible pass/fail evidence for an isolated capture grid."""
    import json

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    cells = tuple(
        measure_cell(manifest_path.parent / cell["path"], cell)
        for cell in payload["cells"]
    )
    expected_count = len(payload["cells"])
    if not cells:
        return _empty_qualification()
    quiet = sum(cell.peak < 1.0e-5 for cell in cells)
    duplicate_pairs = _cross_articulation_duplicates(cells)
    repeat_diversity = _repeat_diversity(cells)
    distances = _articulation_distances(cells)
    velocity_correlations = _velocity_correlations(cells)
    maximum_onset_error = max(abs(cell.onset_error_ms) for cell in cells)
    minimum_distance = min(distances.values()) if distances else 0.0
    minimum_correlation = min(velocity_correlations.values(), default=0.0)
    minimum_unique_repeats = min(repeat_diversity.values(), default=0)
    checks = {
        "cell_count": len(cells) == expected_count and expected_count > 0,
        "no_quiet_cells": quiet == 0,
        "onset_alignment": maximum_onset_error <= 25.0,
        "no_cross_articulation_duplicates": duplicate_pairs == 0,
        "repeat_diversity": minimum_unique_repeats >= 2,
        "articulation_separability": minimum_distance >= 0.5,
        "velocity_ordering": minimum_correlation >= 0.65,
    }
    return {
        "schema": 1,
        "accepted": all(checks.values()),
        "checks": checks,
        "summary": {
            "cell_count": len(cells),
            "quiet_cell_count": quiet,
            "maximum_onset_error_ms": maximum_onset_error,
            "cross_articulation_duplicate_pairs": duplicate_pairs,
            "minimum_unique_repeats_per_condition": minimum_unique_repeats,
            "minimum_articulation_distance_db": minimum_distance,
            "minimum_velocity_spearman": minimum_correlation,
        },
        "articulation_distances_db": distances,
        "velocity_spearman": velocity_correlations,
        "unique_repeats_per_condition": repeat_diversity,
        "cells": [asdict(cell) for cell in cells],
    }


def _empty_qualification() -> dict:
    checks = {
        "cell_count": False,
        "no_quiet_cells": False,
        "onset_alignment": False,
        "no_cross_articulation_duplicates": False,
        "repeat_diversity": False,
        "articulation_separability": False,
        "velocity_ordering": False,
    }
    return {
        "schema": 1,
        "accepted": False,
        "checks": checks,
        "summary": {
            "cell_count": 0,
            "quiet_cell_count": 0,
            "maximum_onset_error_ms": None,
            "cross_articulation_duplicate_pairs": 0,
            "minimum_unique_repeats_per_condition": 0,
            "minimum_articulation_distance_db": 0.0,
            "minimum_velocity_spearman": 0.0,
        },
        "articulation_distances_db": {},
        "velocity_spearman": {},
        "unique_repeats_per_condition": {},
        "cells": [],
    }


def _fixed_window(samples: np.ndarray, start: int, count: int) -> np.ndarray:
    result = np.zeros(count, dtype=np.float64)
    available = max(0, min(count, samples.size - start))
    result[:available] = samples[start : start + available]
    return result


def _band_profile(frequencies: np.ndarray, power: np.ndarray) -> np.ndarray:
    edges = np.geomspace(180.0, 20_000.0, 13)
    energies = np.asarray(
        [
            np.sum(power[(frequencies >= low) & (frequencies < high)])
            for low, high in zip(edges[:-1], edges[1:])
        ]
    )
    fractions = energies / max(float(np.sum(energies)), np.finfo(float).tiny)
    return 10.0 * np.log10(np.maximum(fractions, 1.0e-12))


def _cross_articulation_duplicates(cells: tuple[CellQuality, ...]) -> int:
    groups: dict[str, set[str]] = {}
    for cell in cells:
        groups.setdefault(cell.content_hash, set()).add(cell.articulation)
    return sum(
        len(articulations) * (len(articulations) - 1) // 2
        for articulations in groups.values()
    )


def _repeat_diversity(cells: tuple[CellQuality, ...]) -> dict[str, int]:
    conditions = sorted({(cell.articulation, cell.velocity) for cell in cells})
    return {
        f"{articulation}__v{velocity:03d}": len(
            {
                cell.content_hash
                for cell in cells
                if cell.articulation == articulation and cell.velocity == velocity
            }
        )
        for articulation, velocity in conditions
    }


def _articulation_distances(cells: tuple[CellQuality, ...]) -> dict[str, float]:
    names = sorted({cell.articulation for cell in cells})
    profiles = {
        name: np.median(
            [cell.band_profile_db for cell in cells if cell.articulation == name],
            axis=0,
        )
        for name in names
    }
    return {
        f"{left}__{right}": float(
            np.sqrt(np.mean(np.square(profiles[left] - profiles[right])))
        )
        for index, left in enumerate(names)
        for right in names[index + 1 :]
    }


def _velocity_correlations(cells: tuple[CellQuality, ...]) -> dict[str, float]:
    result = {}
    for articulation in sorted({cell.articulation for cell in cells}):
        selected = [cell for cell in cells if cell.articulation == articulation]
        velocities = sorted({cell.velocity for cell in selected})
        levels = [
            np.median(
                [cell.rms_100ms for cell in selected if cell.velocity == velocity]
            )
            for velocity in velocities
        ]
        result[articulation] = _spearman(
            np.asarray(velocities), np.log(np.maximum(levels, 1.0e-12))
        )
    return result


def _spearman(left: np.ndarray, right: np.ndarray) -> float:
    if left.size < 2 or right.size != left.size:
        return 0.0
    left_ranks = np.argsort(np.argsort(left)).astype(np.float64)
    right_ranks = np.argsort(np.argsort(right)).astype(np.float64)
    correlation = float(np.corrcoef(left_ranks, right_ranks)[0, 1])
    return correlation if np.isfinite(correlation) else 0.0
