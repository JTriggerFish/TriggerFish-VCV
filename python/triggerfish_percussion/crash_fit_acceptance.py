"""Named crash calibration gates that cannot be waived by aggregate loss."""

from __future__ import annotations

import numpy as np

from .crash_fit_common import (
    CrashFitCell,
    body_rms,
    fixed_rate,
    render_cell,
    render_sparse,
)
from .crash_fit_features import REGIONS, maximum_late_regrowth_db, perceptual_features
from .crash_fit_modes import (
    modal_passes,
    modal_tail,
    reference_modal_residual,
    stable_modes,
)
from .crash_fit_prefix import causal_audio_quality
from .crash_fit_t60 import (
    MAXIMUM_BAND_T60_ABSOLUTE_LOG,
    MAXIMUM_BAND_T60_RMSE_LOG,
    MINIMUM_COMMON_T60_BANDS,
    t60_diagnostics,
)
from .crash_model import CrashFit
from .distances import modal_distance
from .modes import match_modes


def acceptance_diagnostics(
    cells: tuple[CrashFitCell, ...], fit: CrashFit
) -> dict[str, object]:
    results = []
    passed = True
    for cell in cells:
        target = perceptual_features(fixed_rate(cell.reference))
        rendered = render_cell(cell, fit, 10.2)
        candidate = perceptual_features(rendered)
        regions = []
        for region, _, _ in REGIONS:
            level_index = target.names.index(f"level/{region}")
            level_delta = 6.0 * (
                candidate.values[level_index] - target.values[level_index]
            )
            band_indices = [
                index
                for index, name in enumerate(target.names)
                if name.startswith(f"erb/{region}/")
            ]
            shape_delta = 12.0 * (
                candidate.values[band_indices] - target.values[band_indices]
            )
            salient = 12.0 * target.values[band_indices] >= -50.0
            selected_delta = shape_delta[salient]
            if selected_delta.size == 0:
                selected_delta = shape_delta
            shape_rmse = float(np.sqrt(np.mean(np.square(selected_delta))))
            maximum_band_delta = float(np.max(np.abs(selected_delta + level_delta)))
            region_passed = bool(
                abs(level_delta) <= 6.0
                and shape_rmse <= 8.0
                and maximum_band_delta <= 12.0
            )
            passed = passed and region_passed
            regions.append(
                {
                    "region": region,
                    "level_delta_db": float(level_delta),
                    "erb_shape_rmse_db": shape_rmse,
                    "maximum_absolute_erb_delta_db": maximum_band_delta,
                    "passed": region_passed,
                }
            )
        results.append({"cell": cell.label, "regions": regions})
        trajectory_indices = [
            index
            for index, name in enumerate(target.names)
            if name.startswith("trajectory/")
        ]
        trajectory_delta = 12.0 * (
            candidate.values[trajectory_indices] - target.values[trajectory_indices]
        )
        salient = 12.0 * target.values[trajectory_indices] >= -70.0
        selected = trajectory_delta[salient]
        if selected.size == 0:
            selected = trajectory_delta
        trajectory_rmse = float(np.sqrt(np.mean(np.square(selected))))
        trajectory_maximum = float(np.max(np.abs(selected)))
        trajectory_passed = trajectory_rmse <= 10.0 and trajectory_maximum <= 24.0
        passed = passed and trajectory_passed
        results[-1]["trajectory"] = {
            "rmse_db": trajectory_rmse,
            "maximum_absolute_delta_db": trajectory_maximum,
            "passed": trajectory_passed,
        }
        regrowth_db = maximum_late_regrowth_db(rendered)
        regrowth_passed = regrowth_db <= 6.0
        passed = passed and regrowth_passed
        results[-1]["stability"] = {
            "maximum_late_regrowth_db": regrowth_db,
            "passed": regrowth_passed,
        }
        texture = causal_audio_quality(
            rendered, cell.reference, 0.100, candidate_onset_sample=0
        )
        texture_passed = bool(
            texture.fine_spectrum_rmse_db <= 4.0
            and texture.centroid_rmse_octaves <= 0.35
            and texture.rolloff_rmse_octaves <= 0.35
            and texture.flatness_rmse_db <= 2.5
            and texture.crest_rmse_db <= 2.5
            and texture.ridge_ratio_absolute_error <= 0.08
        )
        passed = passed and texture_passed
        results[-1]["initial_texture"] = {
            "fine_spectrum_rmse_db": texture.fine_spectrum_rmse_db,
            "centroid_rmse_octaves": texture.centroid_rmse_octaves,
            "rolloff_rmse_octaves": texture.rolloff_rmse_octaves,
            "flatness_rmse_db": texture.flatness_rmse_db,
            "crest_rmse_db": texture.crest_rmse_db,
            "ridge_ratio_absolute_error": texture.ridge_ratio_absolute_error,
            "passed": texture_passed,
        }
        decay = t60_diagnostics(fixed_rate(cell.reference), rendered)
        passed = passed and bool(decay["passed"])
        results[-1]["body_t60"] = decay
    persistent_modes = _persistent_mode_diagnostics(cells, fit)
    passed = passed and bool(persistent_modes["passed"])
    return {
        "passed": passed,
        "gates": {
            "absolute_region_level_db": 6.0,
            "erb_shape_rmse_db": 8.0,
            "maximum_absolute_erb_delta_db": 12.0,
            "reference_relative_band_floor_db": -50.0,
            "trajectory_rmse_db": 10.0,
            "maximum_absolute_trajectory_delta_db": 24.0,
            "maximum_late_regrowth_db": 6.0,
            "maximum_fine_spectrum_rmse_db": 4.0,
            "maximum_centroid_rmse_octaves": 0.35,
            "maximum_rolloff_rmse_octaves": 0.35,
            "maximum_flatness_rmse_db": 2.5,
            "maximum_crest_rmse_db": 2.5,
            "maximum_ridge_ratio_absolute_error": 0.08,
            "maximum_band_t60_rmse_log": MAXIMUM_BAND_T60_RMSE_LOG,
            "maximum_band_t60_absolute_log": MAXIMUM_BAND_T60_ABSOLUTE_LOG,
            "minimum_common_t60_bands": MINIMUM_COMMON_T60_BANDS,
        },
        "cells": results,
        "persistent_modes": persistent_modes,
    }


def _persistent_mode_diagnostics(
    cells: tuple[CrashFitCell, ...], fit: CrashFit
) -> dict[str, object]:
    passes = modal_passes()
    references = tuple(modal_tail(cell.reference) for cell in cells)
    candidates = tuple(modal_tail(render_sparse(cell, fit, 1.2)) for cell in cells)
    required_hits = max(1, int(np.ceil(0.75 * len(cells))))
    reference = stable_modes(references, 48000, passes, required_hits)
    candidate = stable_modes(candidates, 48000, passes, required_hits)
    matching = match_modes(
        reference,
        candidate,
        frequency_scale_cents=20.0,
        decay_scale_log=0.6,
        amplitude_scale_db=18.0,
        maximum_frequency_error_cents=100.0,
        maximum_cost=8.0,
    )
    distances = {term.name: term.value for term in modal_distance(matching)}
    coverage = len(matching.matches) / max(1, len(reference))
    frequency_rms = float(distances["mode_frequency"])
    presentation_deltas = []
    for cell in cells:
        target = body_rms(reference_modal_residual(cell, fit)[0])
        candidate_level = body_rms(render_sparse(cell, fit, 1.6))
        if target > 0.0 and candidate_level > 0.0:
            presentation_deltas.append(20.0 * np.log10(candidate_level / target))
    presentation_level_rms = float(
        np.sqrt(np.mean(np.square(presentation_deltas)))
        if presentation_deltas
        else np.inf
    )
    passed = bool(
        len(reference) >= 3
        and coverage >= 0.75
        and np.isfinite(frequency_rms)
        and frequency_rms <= 50.0
        and float(distances["mode_decay"]) <= 0.5
        and presentation_level_rms <= 6.0
    )
    return {
        "passed": passed,
        "reference_count": len(reference),
        "candidate_count": len(candidate),
        "matched_count": len(matching.matches),
        "reference_coverage": coverage,
        "frequency_rmse_cents": frequency_rms,
        "decay_rmse_log_seconds": float(distances["mode_decay"]),
        "presentation_level_rmse_db": presentation_level_rms,
        "estimator_level_rmse_db": float(distances["mode_level"]),
        "unmatched_count": int(distances["mode_count"]),
        "reference_frequencies_hz": [mode.frequency_hz for mode in reference],
        "candidate_frequencies_hz": [mode.frequency_hz for mode in candidate],
    }
