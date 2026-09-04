"""Regularized object-profile refinement for one crash-cymbal recording."""

from __future__ import annotations

from dataclasses import asdict, replace
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

from .crash_fit_common import CrashFitCell
from .crash_fit_spectral_objective import (
    PREFIXES,
    PROTECTED_LOSS_TOLERANCE,
    make_spectral_targets,
    prefix_losses,
    prefix_qualities,
    profile_residual,
    quality_passes,
    replace_temporal_parameters,
    temporal_bounds,
    temporal_residual,
    temporal_spectral_parameter_names,
)
from .crash_model import CrashFit


class _CountedResidual:
    def __init__(self, function: Callable[[np.ndarray], np.ndarray]):
        self.function = function
        self.calls = 0

    def __call__(self, values: np.ndarray) -> np.ndarray:
        self.calls += 1
        return self.function(values)


def refine_initial_spectral_profile(
    cell: CrashFitCell,
    initial: CrashFit,
    maximum_evaluations: int = 800,
    progress: Callable[[str], None] | None = None,
) -> tuple[CrashFit, dict[str, object]]:
    """Fit one 100 ms private-corpus anchor without moving source balance."""
    diagnostic_evaluations = 8
    if maximum_evaluations < 128 + diagnostic_evaluations:
        raise ValueError("spectral-profile refinement needs at least 136 renders")
    targets = make_spectral_targets(cell)
    optimizer_budget = maximum_evaluations - diagnostic_evaluations
    profile_budget = optimizer_budget * 11 // 20
    profile, profile_stage = _fit_profile(
        cell, initial, targets, profile_budget, progress
    )
    temporal, temporal_stage = _fit_temporal(
        cell, profile, targets, optimizer_budget - profile_budget, progress
    )
    quality = prefix_qualities(cell, temporal)
    accepted = _candidate_passes(quality, profile_stage, temporal_stage)
    profile_stage["accepted"] = _protected_loss_passes(profile_stage)
    temporal_stage["accepted"] = accepted
    temporal_stage["candidate_absolute_quality"] = _quality_summary(cell, quality)
    fitted = temporal if accepted else initial
    return fitted, _fit_diagnostics(
        fitted,
        cell,
        targets,
        maximum_evaluations,
        profile_stage,
        temporal_stage,
        accepted,
        diagnostic_evaluations,
    )


def _fit_profile(cell, initial, targets, budget, progress):
    if progress:
        progress("regularized modal-energy profile refinement")
    start = np.asarray(initial.sparse_amplitude, dtype=np.float64)

    def make_fit(values):
        return replace(initial, sparse_amplitude=tuple(values.tolist()))

    residual = _CountedResidual(
        lambda values: profile_residual(cell, make_fit(values), targets, start)
    )
    result = _solve(residual, start, 0.0, 8.0, budget)
    candidate = make_fit(result.x)
    return candidate, _stage_diagnostics(
        "initial-modal-energy-profile",
        cell,
        initial,
        candidate,
        targets,
        residual.calls,
        result.cost,
        required=False,
    )


def _fit_temporal(cell, initial, targets, budget, progress):
    if progress:
        progress("fixed-balance temporal spectrum refinement")
    names, lower, upper, start = temporal_bounds(initial)

    def make_fit(values):
        return replace_temporal_parameters(initial, names, values)

    residual = _CountedResidual(
        lambda values: temporal_residual(
            cell, make_fit(values), targets, values, start, lower, upper
        )
    )
    result = _solve(residual, start, lower, upper, budget)
    candidate = make_fit(result.x)
    return candidate, _stage_diagnostics(
        "initial-temporal-spectrum",
        cell,
        initial,
        candidate,
        targets,
        residual.calls,
        result.cost,
        required=True,
    )


def _solve(residual, start, lower, upper, budget):
    evaluations_per_step = start.size + 1
    steps = max(4, budget // evaluations_per_step)
    return least_squares(
        residual,
        start,
        bounds=(lower, upper),
        max_nfev=steps,
        x_scale="jac",
        diff_step=0.015,
        ftol=1.0e-5,
        xtol=1.0e-5,
        gtol=1.0e-5,
    )


def _candidate_passes(quality, profile_stage, temporal_stage):
    protected = all(quality_passes(quality[prefix], prefix) for prefix in PREFIXES[:-1])
    return bool(
        protected
        and _protected_loss_passes(profile_stage)
        and _protected_loss_passes(temporal_stage)
        and quality_passes(quality[0.1], 0.1)
    )


def _protected_loss_passes(stage):
    return stage["worst_protected_loss_ratio"] <= 1.0 + PROTECTED_LOSS_TOLERANCE


def _fit_diagnostics(
    fitted, cell, targets, budget, profile, temporal, accepted, diagnostic_evaluations
):
    losses = prefix_losses(cell, fitted, targets)
    return {
        "method": "causal-initial-spectrum-least-squares-v1",
        "policy": "fixed-balance-object-profile",
        "policy_configuration": {
            "protected_acceptance_tolerance": PROTECTED_LOSS_TOLERANCE
        },
        "prefix_seconds": list(PREFIXES),
        "maximum_evaluations": budget,
        "actual_render_evaluations": (
            profile["optimizer_evaluations"]
            + temporal["optimizer_evaluations"]
            + diagnostic_evaluations
        ),
        "requested_final_prefix_seconds": 0.1,
        "start_stage": "initial-modal-energy-profile",
        "seed_prefix_seconds": [0.004, 0.015],
        "workers": 1,
        "completed": accepted,
        "blocked_stage": None if accepted else "initial-temporal-spectrum",
        "stages": [profile, temporal],
        "final_prefix_losses": _string_keys(losses),
    }


def _stage_diagnostics(name, cell, baseline, candidate, targets, calls, cost, required):
    baseline_losses = prefix_losses(cell, baseline, targets)
    candidate_losses = prefix_losses(cell, candidate, targets)
    quality = prefix_qualities(cell, candidate)
    return {
        "stage": name,
        "end_seconds": 0.1,
        "evaluations": calls,
        "optimizer_evaluations": calls,
        "least_squares_cost": float(cost),
        "baseline_prefix_losses": _string_keys(baseline_losses),
        "candidate_prefix_losses": _string_keys(candidate_losses),
        "worst_protected_loss_ratio": max(
            candidate_losses[prefix] / baseline_losses[prefix]
            for prefix in PREFIXES[:-1]
        ),
        "worst_current_loss_ratio": candidate_losses[0.1] / baseline_losses[0.1],
        "accepted": True,
        "candidate_absolute_quality": _quality_summary(cell, quality, required),
        "candidate_cumulative_absolute_quality": [
            {"prefix_seconds": prefix, "passed": quality_passes(item, prefix)}
            for prefix, item in quality.items()
        ],
    }


def _quality_summary(cell, qualities, required=True):
    passed = all(quality_passes(item, prefix) for prefix, item in qualities.items())
    final = qualities[0.1]
    values = asdict(final)
    values.update({"cell": cell.label, "passed": quality_passes(final, 0.1)})
    return {"passed": passed, "required": required, "cells": [values]}


def _string_keys(values):
    return {f"{prefix:.3f}": value for prefix, value in values.items()}
