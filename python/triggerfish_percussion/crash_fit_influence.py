"""Finite-difference checks of when crash parameters affect rendered audio."""

from __future__ import annotations

import numpy as np

from .crash_fit_common import CrashFitCell, render_cell
from .crash_fit_parameters import (
    ALL_CAUSAL_PARAMETERS,
    FitParameter,
    fit_parameter_value,
    replace_fit_parameters,
)
from .crash_fit_prefix import CAUSAL_PREFIX_SECONDS
from .crash_model import CrashFit


def parameter_influence_diagnostics(
    cell: CrashFitCell,
    fit: CrashFit,
    parameters: tuple[FitParameter, ...] = ALL_CAUSAL_PARAMETERS,
    seconds: float = 8.0,
) -> dict[str, object]:
    """Measure finite-difference energy in disjoint causal time regions."""
    duration = max(seconds, CAUSAL_PREFIX_SECONDS[0])
    baseline = render_cell(cell, fit, duration).samples
    results = []
    for specification in parameters:
        perturbed_value, delta, range_fraction = _perturbation(fit, specification)
        candidate = render_cell(
            cell,
            replace_fit_parameters(fit, {specification.name: perturbed_value}),
            duration,
        ).samples
        # Normalize by the fraction of the declared parameter range, not by
        # raw units. Hertz, seconds and dimensionless gains are then comparable.
        derivative = (candidate - baseline) / range_fraction
        regions = _region_levels(derivative, baseline, 48000, duration)
        detectable = [
            region
            for region in regions
            if region["range_normalized_relative_db"] > -80.0
        ]
        results.append(
            {
                "parameter": specification.name,
                "delta": delta,
                "range_fraction": range_fraction,
                "earliest_detectable_seconds": (
                    detectable[0]["start_seconds"] if detectable else None
                ),
                "regions": regions,
            }
        )
    return {
        "cell": cell.label,
        "strength": cell.strength,
        "location": cell.location,
        "hardness": cell.hardness,
        "detection_floor_db": -80.0,
        "parameters": results,
    }


def _perturbation(
    fit: CrashFit, specification: FitParameter
) -> tuple[float, float, float]:
    value = fit_parameter_value(fit, specification.name)
    span = specification.upper - specification.lower
    step = max(specification.relative_step * span, 1.0e-4 * span)
    perturbed = min(specification.upper, value + step)
    if perturbed == value:
        perturbed = max(specification.lower, value - step)
    delta = perturbed - value
    return perturbed, delta, abs(delta) / span


def _region_levels(
    derivative: np.ndarray,
    baseline: np.ndarray,
    sample_rate: int,
    seconds: float,
) -> list[dict[str, float]]:
    ends = [value for value in CAUSAL_PREFIX_SECONDS if value < seconds]
    ends.append(seconds)
    start = 0.0
    result = []
    for end in ends:
        first = min(derivative.size, round(start * sample_rate))
        last = min(derivative.size, round(end * sample_rate))
        difference_energy = float(np.sum(np.square(derivative[first:last])))
        baseline_energy = float(np.sum(np.square(baseline[first:last])))
        relative = 10.0 * np.log10(
            max(difference_energy, 1.0e-30) / max(baseline_energy, 1.0e-30)
        )
        result.append(
            {
                "start_seconds": start,
                "end_seconds": end,
                "range_normalized_relative_db": relative,
            }
        )
        start = end
    return result
