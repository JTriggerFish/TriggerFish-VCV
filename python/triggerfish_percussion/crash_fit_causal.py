"""Cumulative crash fitting with explicit prefix trade-off policies."""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from time import perf_counter

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc

from .crash_fit_common import CrashFitCell, render_cell
from .crash_fit_parameters import (
    CAUSAL_STAGES,
    CausalStage,
    fit_parameter_value,
    replace_fit_parameters,
)
from .crash_fit_prefix import (
    CAUSAL_PREFIX_SECONDS,
    CausalFeatures,
    PrefixQuality,
    causal_feature_loss,
    causal_prefix_quality,
    causal_prefix_features,
)
from .crash_fit_texture import (
    TextureDescriptor,
    compare_texture_descriptors,
    texture_descriptor,
)
from .crash_model import CrashFit


@dataclass(frozen=True)
class CausalFitPolicy:
    """Control how a longer prefix may trade accuracy against earlier prefixes."""

    name: str
    prefix_weights: tuple[tuple[float, float], ...] = ()
    protected_objective_tolerance: float = 0.02
    protected_acceptance_tolerance: float = 0.025
    protected_penalty: float = 1000.0
    protected_quality_barrier: float = 1.0e6
    absolute_penalty: float = 1000.0
    require_absolute_gate: bool = True
    stop_on_absolute_failure: bool = True

    def weight(self, prefix: float) -> float:
        return dict(self.prefix_weights).get(prefix, 1.0)


STRICT_CAUSAL_POLICY = CausalFitPolicy("strict")

# The 100 ms feature already includes the attack. These smaller weights prevent
# its 0-4 and 0-15 ms content from being counted three times, while the loose
# per-cell guard still rejects catastrophic attack regressions.
FIRST_100MS_TRADEOFF_POLICY = CausalFitPolicy(
    "first-100ms-tradeoff",
    prefix_weights=((0.004, 0.05), (0.015, 0.15), (0.100, 1.0)),
    protected_objective_tolerance=0.35,
    protected_acceptance_tolerance=0.35,
    protected_penalty=20.0,
    absolute_penalty=1000.0,
    require_absolute_gate=True,
    stop_on_absolute_failure=True,
)

CAUSAL_FIT_POLICIES = {
    policy.name: policy
    for policy in (STRICT_CAUSAL_POLICY, FIRST_100MS_TRADEOFF_POLICY)
}


def fit_causal_model(
    cells: tuple[CrashFitCell, ...],
    initial: CrashFit,
    maximum_evaluations: int = 800,
    seed: int = 1234,
    stages: tuple[CausalStage, ...] = CAUSAL_STAGES,
    start_stage: str | None = None,
    policy: CausalFitPolicy = STRICT_CAUSAL_POLICY,
    progress: Callable[[str], None] | None = None,
    workers: int = 1,
) -> tuple[CrashFit, dict[str, object]]:
    """Fit progressively longer prefixes without invalidating earlier ones."""
    if not cells:
        raise ValueError("causal crash fitting needs at least one cell")
    if not stages:
        raise ValueError("causal crash fitting needs at least one stage")
    if workers < 1:
        raise ValueError("causal crash fitting needs at least one worker")
    prefixes = tuple(sorted({stage.end_seconds for stage in stages}))
    targets = _target_features(cells, prefixes)
    texture_targets = _target_texture_descriptors(cells, prefixes)
    start_index = _start_index(stages, start_stage)
    active_stages = stages[start_index:]
    budgets = list(_stage_budgets(active_stages, maximum_evaluations))
    current = initial
    prefix_limits: dict[tuple[int, float], float] = {}
    stage_diagnostics = []
    frozen_prefixes = tuple(
        sorted(
            {
                stage.end_seconds
                for stage in stages[:start_index]
                if stage.requires_acceptance_gate
            }
        )
    )
    attempted_prefixes: set[float] = set(frozen_prefixes)
    if frozen_prefixes:
        frozen_losses = _cell_prefix_losses(cells, current, frozen_prefixes, targets)
        prefix_limits.update(frozen_losses)
        _validate_frozen_stages(
            cells,
            current,
            stages[:start_index],
            targets,
            texture_targets,
            start_stage,
        )
    blocked_stage = None
    for local_index, (stage, budget) in enumerate(zip(active_stages, budgets)):
        index = start_index + local_index
        if budget < 8:
            blocked_stage = f"{stage.name}: evaluation budget exhausted"
            break
        if progress:
            progress(
                f"causal stage {index + 1}/{len(stages)}: {stage.name} "
                f"(0-{stage.end_seconds:g} s, {budget} evaluations)"
            )
        stage_start = perf_counter()
        cumulative = tuple(prefix for prefix in prefixes if prefix <= stage.end_seconds)
        current, diagnostics, accepted_losses = _fit_stage(
            cells,
            targets,
            texture_targets,
            current,
            stage,
            cumulative,
            prefix_limits,
            budget,
            seed + 7919 * index,
            policy,
            workers,
        )
        for key, loss in accepted_losses.items():
            prefix_limits[key] = min(prefix_limits.get(key, loss), loss)
        attempted_prefixes.add(stage.end_seconds)
        diagnostics["elapsed_seconds"] = perf_counter() - stage_start
        stage_diagnostics.append(diagnostics)
        if progress:
            status = "accepted" if diagnostics["accepted"] else "rejected"
            progress(
                f"causal stage {stage.name}: {status}; "
                f"{diagnostics['evaluations']} evaluations in "
                f"{diagnostics['elapsed_seconds']:.1f} s"
            )
        if (
            policy.stop_on_absolute_failure
            and stage.requires_quality
            and stage.requires_acceptance_gate
            and not diagnostics["absolute_quality"]["passed"]
        ):
            retry_budget = min(budget, sum(budgets[local_index + 1 :]))
            if retry_budget >= 8:
                _consume_later_budgets(budgets, local_index, retry_budget)
                if progress:
                    progress(
                        f"causal stage {stage.name}: absolute gate failed; "
                        f"reallocating {retry_budget} later-stage evaluations"
                    )
                retry_start = perf_counter()
                retry_stage = _cumulative_retry_stage(stages[: index + 1], stage)
                current, retry, accepted_losses = _fit_stage(
                    cells,
                    targets,
                    texture_targets,
                    current,
                    retry_stage,
                    cumulative,
                    prefix_limits,
                    retry_budget,
                    seed + 104729 * (index + 1),
                    policy,
                    workers,
                )
                for key, loss in accepted_losses.items():
                    prefix_limits[key] = min(prefix_limits.get(key, loss), loss)
                retry["retry"] = True
                retry["stage"] = stage.name
                retry["reopened_parameters"] = [
                    item.name for item in retry_stage.parameters
                ]
                retry["elapsed_seconds"] = perf_counter() - retry_start
                stage_diagnostics.append(retry)
                diagnostics = retry
                if progress:
                    quality = (
                        "passed" if retry["absolute_quality"]["passed"] else "failed"
                    )
                    progress(
                        f"causal stage {stage.name} retry: absolute gate {quality} "
                        f"in {retry['elapsed_seconds']:.1f} s"
                    )
            if (
                policy.stop_on_absolute_failure
                and not diagnostics["absolute_quality"]["passed"]
            ):
                blocked_stage = stage.name
                if progress:
                    progress(
                        f"causal fit stopped at {stage.name}: active prefix is not matched"
                    )
                break
    fitted_prefixes = tuple(sorted(attempted_prefixes))
    losses = causal_prefix_losses(cells, current, fitted_prefixes, targets)
    return current, {
        "method": "cumulative-causal-prefix-v3",
        "policy": policy.name,
        "policy_configuration": {
            "prefix_weights": {
                f"{prefix:.3f}": weight for prefix, weight in policy.prefix_weights
            },
            "protected_objective_tolerance": policy.protected_objective_tolerance,
            "protected_acceptance_tolerance": policy.protected_acceptance_tolerance,
            "protected_quality_barrier": policy.protected_quality_barrier,
            "absolute_gate_required": policy.require_absolute_gate,
        },
        "prefix_seconds": list(fitted_prefixes),
        "maximum_evaluations": maximum_evaluations,
        "requested_final_prefix_seconds": max(stage.end_seconds for stage in stages),
        "start_stage": start_stage,
        "seed_prefix_seconds": list(frozen_prefixes),
        "workers": workers,
        "completed": blocked_stage is None,
        "blocked_stage": blocked_stage,
        "stages": stage_diagnostics,
        "final_prefix_losses": _string_keys(losses),
    }


def causal_prefix_losses(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    prefixes: tuple[float, ...] = CAUSAL_PREFIX_SECONDS,
    targets: dict[tuple[int, float], CausalFeatures] | None = None,
) -> dict[float, float]:
    target_features = targets or _target_features(cells, prefixes)
    cell_losses = _cell_prefix_losses(cells, fit, prefixes, target_features)
    return {
        prefix: float(
            np.mean(
                [cell_losses[(cell_index, prefix)] for cell_index in range(len(cells))]
            )
        )
        for prefix in prefixes
    }


def _cell_prefix_losses(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    prefixes: tuple[float, ...],
    targets: dict[tuple[int, float], CausalFeatures],
) -> dict[tuple[int, float], float]:
    return _cell_prefix_evaluation(cells, fit, prefixes, targets)[0]


def _cell_prefix_evaluation(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    prefixes: tuple[float, ...],
    targets: dict[tuple[int, float], CausalFeatures],
    quality_prefixes: tuple[float, ...] = (),
    texture_targets: dict[tuple[int, float], TextureDescriptor] | None = None,
) -> tuple[dict[tuple[int, float], float], dict[tuple[int, float], PrefixQuality]]:
    losses = {}
    qualities = {}
    duration = max(prefixes) + 0.05
    for cell_index, cell in enumerate(cells):
        rendered = render_cell(cell, fit, duration)
        for prefix in prefixes:
            candidate = causal_prefix_features(rendered, prefix, onset_sample=0)
            target = targets[(cell_index, prefix)]
            losses[(cell_index, prefix)] = causal_feature_loss(candidate, target)
            if prefix in quality_prefixes:
                base = causal_prefix_quality(candidate, target)
                if prefix < 0.015:
                    qualities[(cell_index, prefix)] = base
                    continue
                target_texture = (texture_targets or {})[(cell_index, prefix)]
                texture = compare_texture_descriptors(
                    texture_descriptor(rendered, prefix, onset_sample=0),
                    target_texture,
                )
                qualities[(cell_index, prefix)] = replace(
                    base,
                    fine_spectrum_rmse_db=texture.fine_spectrum_rmse_db,
                    centroid_rmse_octaves=texture.centroid_rmse_octaves,
                    rolloff_rmse_octaves=texture.rolloff_rmse_octaves,
                    flatness_rmse_db=texture.flatness_rmse_db,
                    crest_rmse_db=texture.crest_rmse_db,
                    ridge_ratio_absolute_error=texture.ridge_ratio_absolute_error,
                )
    return losses, qualities


def prefix_gate_passes(
    baseline: dict[float, float],
    candidate: dict[float, float],
    protected_prefixes: tuple[float, ...],
    relative_tolerance: float = 0.02,
) -> bool:
    return all(
        candidate[prefix] <= baseline[prefix] * (1.0 + relative_tolerance) + 1.0e-9
        for prefix in protected_prefixes
    )


def _fit_stage(
    cells: tuple[CrashFitCell, ...],
    targets: dict[tuple[int, float], CausalFeatures],
    texture_targets: dict[tuple[int, float], TextureDescriptor],
    initial: CrashFit,
    stage: CausalStage,
    prefixes: tuple[float, ...],
    prefix_limits: dict[tuple[int, float], float],
    maximum_evaluations: int,
    seed: int,
    policy: CausalFitPolicy,
    workers: int,
) -> tuple[CrashFit, dict[str, object], dict[tuple[int, float], float]]:
    protected = tuple(prefix for prefix in prefixes if prefix < stage.end_seconds)
    baseline_cells = _cell_prefix_losses(cells, initial, prefixes, targets)
    baseline = _mean_prefix_losses(baseline_cells, len(cells), prefixes)
    names = tuple(item.name for item in stage.parameters)
    bounds = np.asarray(
        [(item.lower, item.upper) for item in stage.parameters], dtype=np.float64
    )
    start = _validated_parameter_start(
        np.asarray(
            [fit_parameter_value(initial, name) for name in names],
            dtype=np.float64,
        ),
        bounds,
        names,
    )
    cache: dict[bytes, tuple[float, dict[tuple[int, float], float]]] = {}

    def evaluate(
        vector: np.ndarray,
    ) -> tuple[float, dict[tuple[int, float], float]]:
        bounded = np.clip(vector, bounds[:, 0], bounds[:, 1])
        key = bounded.astype(np.float64).tobytes()
        if key in cache:
            return cache[key]
        fit = replace_fit_parameters(initial, dict(zip(names, bounded.tolist())))
        losses, qualities = _cell_prefix_evaluation(
            cells,
            fit,
            prefixes,
            targets,
            prefixes if stage.requires_quality else (),
            texture_targets,
        )
        normalized = [
            policy.weight(key[1]) * losses[key] / max(baseline_cells[key], 1.0e-8)
            for key in baseline_cells
        ]
        regressions = [
            max(
                0.0,
                losses[(cell_index, prefix)]
                / max(
                    prefix_limits.get(
                        (cell_index, prefix), baseline_cells[(cell_index, prefix)]
                    ),
                    1.0e-8,
                )
                - (1.0 + policy.protected_objective_tolerance),
            )
            for prefix in protected
            for cell_index in range(len(cells))
        ]
        violations = (
            _quality_violations(qualities, stage) if stage.requires_quality else ()
        )
        protected_quality_passed = _protected_quality_passes(
            qualities, stage, protected
        )
        score = float(
            np.mean(normalized)
            + policy.protected_penalty * np.sum(np.square(regressions))
            + (0.0 if protected_quality_passed else policy.protected_quality_barrier)
            + policy.absolute_penalty * np.sum(np.square(violations))
        )
        cache[key] = score, losses
        return score, losses

    def objective(vector: np.ndarray) -> float:
        return evaluate(vector)[0]

    population = _initial_population(start, bounds, maximum_evaluations, seed)
    feasibility_probes = 0
    if protected and stage.requires_quality:
        population, feasibility_probes = _feasible_local_population(
            start,
            bounds,
            population.shape[0],
            objective,
            policy.protected_quality_barrier,
            seed + 17,
        )
    reserved_local = max(4 * start.size, maximum_evaluations // 5)
    if feasibility_probes:
        remaining_global = max(0, maximum_evaluations - reserved_local - len(cache))
        generations = remaining_global // population.shape[0]
    else:
        global_budget = max(population.shape[0], maximum_evaluations - reserved_local)
        generations = max(0, global_budget // population.shape[0] - 1)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        global_fit = differential_evolution(
            objective,
            bounds,
            init=population,
            x0=start,
            seed=seed,
            maxiter=generations,
            polish=False,
            tol=0.0,
            atol=0.0,
            updating="deferred" if workers > 1 else "immediate",
            workers=executor.map if workers > 1 else 1,
        )
    global_evaluations = len(cache)
    remaining = max(0, maximum_evaluations - len(cache))
    result = global_fit.x
    global_score = objective(result)
    local_score = global_score
    if remaining >= 4:
        polished = minimize(
            objective,
            result,
            method="Powell",
            bounds=bounds,
            options={"maxfev": remaining, "xtol": 1.0e-3, "ftol": 1.0e-4},
        )
        local_score = objective(polished.x)
        if local_score < global_score:
            result = polished.x
    local_evaluations = len(cache) - global_evaluations
    candidate = replace_fit_parameters(initial, dict(zip(names, result.tolist())))
    candidate_cells, candidate_qualities = _cell_prefix_evaluation(
        cells,
        candidate,
        prefixes,
        targets,
        prefixes if stage.requires_quality else (),
        texture_targets,
    )
    candidate_losses = _mean_prefix_losses(candidate_cells, len(cells), prefixes)
    baseline_objective = evaluate(start)[0]
    candidate_objective = evaluate(result)[0]
    protected_limits = {
        (cell_index, prefix): prefix_limits.get(
            (cell_index, prefix), baseline_cells[(cell_index, prefix)]
        )
        for prefix in protected
        for cell_index in range(len(cells))
    }
    current_limits = {
        (cell_index, stage.end_seconds): baseline_cells[(cell_index, stage.end_seconds)]
        for cell_index in range(len(cells))
    }
    gate_passed = _cell_gate_passes(
        protected_limits,
        candidate_cells,
        relative_tolerance=policy.protected_acceptance_tolerance,
    )
    protected_quality_passed = _protected_quality_passes(
        candidate_qualities, stage, protected
    )
    cumulative_quality_passed = bool(
        not stage.requires_quality
        or not any(_quality_violations(candidate_qualities, stage))
    )
    cumulative_quality = _cumulative_quality_diagnostics(
        cells, prefixes, candidate_qualities, stage
    )
    improved = candidate_objective < baseline_objective
    quality_accepted = bool(
        cumulative_quality_passed
        or not policy.require_absolute_gate
        or not stage.requires_acceptance_gate
    )
    accepted = bool(
        gate_passed and protected_quality_passed and quality_accepted and improved
    )
    fitted = candidate if accepted else initial
    accepted_losses = candidate_cells if accepted else baseline_cells
    candidate_quality = _quality_diagnostics(
        cells, candidate, stage, targets, texture_targets
    )
    quality = (
        candidate_quality
        if accepted
        else _quality_diagnostics(cells, initial, stage, targets, texture_targets)
    )
    return (
        fitted,
        {
            "stage": stage.name,
            "policy": policy.name,
            "protected_acceptance_tolerance": (policy.protected_acceptance_tolerance),
            "end_seconds": stage.end_seconds,
            "parameters": list(names),
            "candidate_parameters": dict(zip(names, result.tolist())),
            "evaluation_budget": maximum_evaluations,
            "evaluations": len(cache),
            "global_evaluations": global_evaluations,
            "local_evaluations": local_evaluations,
            "feasible_population_size": population.shape[0],
            "feasibility_probes": feasibility_probes,
            "global_objective": global_score,
            "local_objective": local_score,
            "local_polish_accepted": local_score < global_score,
            "baseline_prefix_losses": _string_keys(baseline),
            "protected_prefix_limits": _mean_prefix_losses_as_strings(
                protected_limits, len(cells), protected
            ),
            "candidate_prefix_losses": _string_keys(candidate_losses),
            "baseline_composite_objective": baseline_objective,
            "candidate_composite_objective": candidate_objective,
            "protected_prefixes": list(protected),
            "prefix_gate_passed": gate_passed,
            "protected_absolute_gate_passed": protected_quality_passed,
            "cumulative_absolute_gate_passed": cumulative_quality_passed,
            "worst_protected_loss_ratio": _worst_ratio(
                protected_limits, candidate_cells
            ),
            "worst_current_loss_ratio": _worst_ratio(current_limits, candidate_cells),
            "current_prefix_improved": (
                candidate_losses[stage.end_seconds] < baseline[stage.end_seconds]
            ),
            "composite_objective_improved": improved,
            "accepted": accepted,
            "candidate_absolute_quality": candidate_quality,
            "candidate_cumulative_absolute_quality": cumulative_quality,
            "absolute_quality": quality,
        },
        accepted_losses,
    )


def _quality_diagnostics(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    stage: CausalStage,
    targets: dict[tuple[int, float], CausalFeatures],
    texture_targets: dict[tuple[int, float], TextureDescriptor],
) -> dict[str, object]:
    duration = stage.end_seconds + 0.05
    results = []
    passed = True
    for cell_index, cell in enumerate(cells):
        rendered = render_cell(cell, fit, duration)
        candidate = causal_prefix_features(rendered, stage.end_seconds, onset_sample=0)
        quality = causal_prefix_quality(
            candidate, targets[(cell_index, stage.end_seconds)]
        )
        if stage.end_seconds >= 0.015:
            texture = compare_texture_descriptors(
                texture_descriptor(rendered, stage.end_seconds, onset_sample=0),
                texture_targets[(cell_index, stage.end_seconds)],
            )
            quality = replace(
                quality,
                fine_spectrum_rmse_db=texture.fine_spectrum_rmse_db,
                centroid_rmse_octaves=texture.centroid_rmse_octaves,
                rolloff_rmse_octaves=texture.rolloff_rmse_octaves,
                flatness_rmse_db=texture.flatness_rmse_db,
                crest_rmse_db=texture.crest_rmse_db,
                ridge_ratio_absolute_error=texture.ridge_ratio_absolute_error,
            )
        cell_passed = _quality_passes(quality, stage, stage.end_seconds)
        passed = passed and cell_passed
        results.append(
            {
                "cell": cell.label,
                "envelope_rmse_db": quality.envelope_rmse_db,
                "envelope_maximum_absolute_db": (quality.envelope_maximum_absolute_db),
                "spectral_rmse_db": quality.spectral_rmse_db,
                "spectral_p95_absolute_db": quality.spectral_p95_absolute_db,
                "fine_spectrum_rmse_db": quality.fine_spectrum_rmse_db,
                "centroid_rmse_octaves": quality.centroid_rmse_octaves,
                "rolloff_rmse_octaves": quality.rolloff_rmse_octaves,
                "flatness_rmse_db": quality.flatness_rmse_db,
                "crest_rmse_db": quality.crest_rmse_db,
                "ridge_ratio_absolute_error": quality.ridge_ratio_absolute_error,
                "passed": cell_passed,
            }
        )
    return {
        "passed": passed,
        "required": stage.requires_quality and stage.requires_acceptance_gate,
        "gates": {
            "maximum_envelope_rmse_db": stage.maximum_envelope_rmse_db,
            "maximum_envelope_absolute_db": stage.maximum_envelope_absolute_db,
            "maximum_spectral_rmse_db": stage.maximum_spectral_rmse_db,
            "maximum_spectral_p95_absolute_db": stage.maximum_spectral_p95_db,
            "maximum_fine_spectrum_rmse_db": stage.maximum_fine_spectrum_rmse_db,
            "maximum_centroid_rmse_octaves": stage.maximum_centroid_rmse_octaves,
            "maximum_rolloff_rmse_octaves": stage.maximum_rolloff_rmse_octaves,
            "maximum_flatness_rmse_db": stage.maximum_flatness_rmse_db,
            "maximum_crest_rmse_db": stage.maximum_crest_rmse_db,
            "maximum_ridge_ratio_error": stage.maximum_ridge_ratio_error,
        },
        "cells": results,
    }


def _quality_violations(
    qualities: dict[tuple[int, float], PrefixQuality], stage: CausalStage
) -> tuple[float, ...]:
    result = []
    for (_, prefix), quality in qualities.items():
        result.extend(
            (
                max(
                    0.0,
                    quality.envelope_rmse_db / stage.maximum_envelope_rmse_db - 1.0,
                ),
                max(
                    0.0,
                    quality.envelope_maximum_absolute_db
                    / stage.maximum_envelope_absolute_db
                    - 1.0,
                ),
                max(
                    0.0,
                    quality.spectral_rmse_db / stage.maximum_spectral_rmse_db - 1.0,
                ),
                max(
                    0.0,
                    quality.spectral_p95_absolute_db / stage.maximum_spectral_p95_db
                    - 1.0,
                ),
            )
        )
        if prefix >= 0.015:
            result.extend(_texture_quality_violations(quality, stage))
    return tuple(result)


def _protected_quality_passes(
    qualities: dict[tuple[int, float], PrefixQuality],
    stage: CausalStage,
    protected_prefixes: tuple[float, ...],
) -> bool:
    """Keep every previously accepted absolute gate during continuation."""
    if not stage.requires_quality:
        return True
    protected = set(protected_prefixes)
    return all(
        _quality_passes(quality, stage, prefix)
        for (_, prefix), quality in qualities.items()
        if prefix in protected
    )


def _cumulative_quality_diagnostics(
    cells: tuple[CrashFitCell, ...],
    prefixes: tuple[float, ...],
    qualities: dict[tuple[int, float], PrefixQuality],
    stage: CausalStage,
) -> list[dict[str, object]]:
    if not stage.requires_quality:
        return []
    result = []
    for prefix in prefixes:
        values = [qualities[(cell_index, prefix)] for cell_index in range(len(cells))]
        envelope = max(item.envelope_rmse_db for item in values)
        envelope_maximum = max(item.envelope_maximum_absolute_db for item in values)
        spectral = max(item.spectral_rmse_db for item in values)
        p95 = max(item.spectral_p95_absolute_db for item in values)
        result.append(
            {
                "prefix_seconds": prefix,
                "maximum_envelope_rmse_db": envelope,
                "maximum_envelope_absolute_db": envelope_maximum,
                "maximum_spectral_rmse_db": spectral,
                "maximum_spectral_p95_absolute_db": p95,
                "maximum_fine_spectrum_rmse_db": max(
                    item.fine_spectrum_rmse_db for item in values
                ),
                "maximum_centroid_rmse_octaves": max(
                    item.centroid_rmse_octaves for item in values
                ),
                "maximum_rolloff_rmse_octaves": max(
                    item.rolloff_rmse_octaves for item in values
                ),
                "maximum_flatness_rmse_db": max(
                    item.flatness_rmse_db for item in values
                ),
                "maximum_crest_rmse_db": max(item.crest_rmse_db for item in values),
                "maximum_ridge_ratio_absolute_error": max(
                    item.ridge_ratio_absolute_error for item in values
                ),
                "passed": all(_quality_passes(item, stage, prefix) for item in values),
            }
        )
    return result


def _quality_passes(
    quality: PrefixQuality, stage: CausalStage, prefix_seconds: float
) -> bool:
    base_passes = bool(
        quality.envelope_rmse_db <= stage.maximum_envelope_rmse_db
        and quality.envelope_maximum_absolute_db <= stage.maximum_envelope_absolute_db
        and quality.spectral_rmse_db <= stage.maximum_spectral_rmse_db
        and quality.spectral_p95_absolute_db <= stage.maximum_spectral_p95_db
    )
    if prefix_seconds < 0.015:
        return base_passes
    return base_passes and not any(_texture_quality_violations(quality, stage))


def _texture_quality_violations(
    quality: PrefixQuality, stage: CausalStage
) -> tuple[float, ...]:
    return (
        max(
            0.0,
            quality.fine_spectrum_rmse_db / stage.maximum_fine_spectrum_rmse_db - 1.0,
        ),
        max(
            0.0,
            quality.centroid_rmse_octaves / stage.maximum_centroid_rmse_octaves - 1.0,
        ),
        max(
            0.0,
            quality.rolloff_rmse_octaves / stage.maximum_rolloff_rmse_octaves - 1.0,
        ),
        max(0.0, quality.flatness_rmse_db / stage.maximum_flatness_rmse_db - 1.0),
        max(0.0, quality.crest_rmse_db / stage.maximum_crest_rmse_db - 1.0),
        max(
            0.0,
            quality.ridge_ratio_absolute_error / stage.maximum_ridge_ratio_error - 1.0,
        ),
    )


def _target_features(
    cells: tuple[CrashFitCell, ...], prefixes: tuple[float, ...]
) -> dict[tuple[int, float], CausalFeatures]:
    return {
        (cell_index, prefix): causal_prefix_features(cell.reference, prefix)
        for cell_index, cell in enumerate(cells)
        for prefix in prefixes
    }


def _target_texture_descriptors(
    cells: tuple[CrashFitCell, ...], prefixes: tuple[float, ...]
) -> dict[tuple[int, float], TextureDescriptor]:
    return {
        (cell_index, prefix): texture_descriptor(cell.reference, prefix)
        for cell_index, cell in enumerate(cells)
        for prefix in prefixes
        if prefix >= 0.015
    }


def _cell_gate_passes(
    limits: dict[tuple[int, float], float],
    candidate: dict[tuple[int, float], float],
    relative_tolerance: float = 0.02,
) -> bool:
    return all(
        candidate[key] <= limit * (1.0 + relative_tolerance) + 1.0e-9
        for key, limit in limits.items()
    )


def _worst_ratio(
    limits: dict[tuple[int, float], float],
    candidate: dict[tuple[int, float], float],
) -> float:
    if not limits:
        return 1.0
    return float(
        max(candidate[key] / max(limit, 1.0e-8) for key, limit in limits.items())
    )


def _mean_prefix_losses(
    values: dict[tuple[int, float], float],
    cell_count: int,
    prefixes: tuple[float, ...],
) -> dict[float, float]:
    return {
        prefix: float(
            np.mean([values[(cell_index, prefix)] for cell_index in range(cell_count)])
        )
        for prefix in prefixes
    }


def _mean_prefix_losses_as_strings(
    values: dict[tuple[int, float], float],
    cell_count: int,
    prefixes: tuple[float, ...],
) -> dict[str, float]:
    if not prefixes:
        return {}
    return _string_keys(_mean_prefix_losses(values, cell_count, prefixes))


def _stage_budgets(
    stages: tuple[CausalStage, ...], maximum_evaluations: int
) -> tuple[int, ...]:
    minimum = 8 * len(stages)
    total = max(minimum, maximum_evaluations)
    weights = np.asarray([np.sqrt(len(stage.parameters)) for stage in stages])
    raw = (total - minimum) * weights / np.sum(weights)
    budgets = 8 + np.floor(raw).astype(int)
    for index in range(total - int(np.sum(budgets))):
        budgets[index % len(budgets)] += 1
    return tuple(int(value) for value in budgets)


def _start_index(stages: tuple[CausalStage, ...], start_stage: str | None) -> int:
    if start_stage is None:
        return 0
    for index, stage in enumerate(stages):
        if stage.name == start_stage:
            return index
    raise ValueError(f"unknown causal start stage: {start_stage}")


def _validate_frozen_stages(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    stages: tuple[CausalStage, ...],
    targets: dict[tuple[int, float], CausalFeatures],
    texture_targets: dict[tuple[int, float], TextureDescriptor],
    start_stage: str | None,
) -> None:
    for stage in stages:
        if stage.requires_quality and stage.requires_acceptance_gate:
            quality = _quality_diagnostics(cells, fit, stage, targets, texture_targets)
            if not quality["passed"]:
                raise ValueError(
                    f"seed fit fails {stage.name}; cannot resume at {start_stage}"
                )


def _consume_later_budgets(budgets: list[int], stage_index: int, amount: int) -> None:
    """Deduct a retry allocation from the furthest unstarted stages."""
    remaining = amount
    for index in range(len(budgets) - 1, stage_index, -1):
        consumed = min(remaining, budgets[index])
        budgets[index] -= consumed
        remaining -= consumed
        if remaining == 0:
            return


def _cumulative_retry_stage(
    active_stages: tuple[CausalStage, ...], current: CausalStage
) -> CausalStage:
    """Reopen every control already capable of affecting the failed prefix."""
    parameters = {
        item.name: item
        for stage in active_stages
        if stage.end_seconds <= current.end_seconds
        for item in stage.parameters
    }
    return replace(current, parameters=tuple(parameters.values()))


def _validated_parameter_start(
    start: np.ndarray, bounds: np.ndarray, names: tuple[str, ...]
) -> np.ndarray:
    """Reject undeclared seed clipping before an optimization begins."""
    outside = (start < bounds[:, 0]) | (start > bounds[:, 1])
    if np.any(outside):
        index = int(np.flatnonzero(outside)[0])
        raise ValueError(
            f"fit seed {names[index]}={start[index]:g} lies outside "
            f"[{bounds[index, 0]:g}, {bounds[index, 1]:g}]"
        )
    return start


def _initial_population(
    start: np.ndarray, bounds: np.ndarray, budget: int, seed: int
) -> np.ndarray:
    maximum_count = max(8, min(8 * start.size, budget // 4))
    count = 1 << int(np.floor(np.log2(maximum_count)))
    span = bounds[:, 1] - bounds[:, 0]
    exponent = int(np.log2(count))
    global_unit = qmc.Sobol(start.size, scramble=True, seed=seed).random_base2(exponent)
    population = qmc.scale(global_unit, bounds[:, 0], bounds[:, 1])
    population[0] = start
    local_count = count // 2 - 1
    local_unit = qmc.Sobol(start.size, scramble=True, seed=seed + 1).random_base2(
        exponent - 1
    )[:local_count]
    population[1 : 1 + local_count] = np.clip(
        start + (2.0 * local_unit - 1.0) * 0.2 * span,
        bounds[:, 0],
        bounds[:, 1],
    )
    return population


def _feasible_local_population(
    start: np.ndarray,
    bounds: np.ndarray,
    target_count: int,
    objective: Callable[[np.ndarray], float],
    barrier: float,
    seed: int,
) -> tuple[np.ndarray, int]:
    """Build a deterministic trust-region population inside hard fit gates."""
    accepted = [start.copy()]
    seen = {start.astype(np.float64).tobytes()}
    probes = 1
    if objective(start) >= barrier:
        raise ValueError("constrained population seed violates a protected gate")
    span = bounds[:, 1] - bounds[:, 0]
    dimension = start.size
    for radius_index, radius in enumerate((0.1, 0.05, 0.025, 0.0125, 0.00625)):
        directions = []
        for axis in range(dimension):
            direction = np.zeros(dimension)
            direction[axis] = 1.0
            directions.extend((direction, -direction))
        unit = qmc.Sobol(dimension, scramble=True, seed=seed + radius_index).random(
            2 * max(target_count, dimension)
        )
        directions.extend(2.0 * unit - 1.0)
        for direction in directions:
            candidate = np.clip(
                start + radius * span * direction, bounds[:, 0], bounds[:, 1]
            )
            key = candidate.astype(np.float64).tobytes()
            if key in seen:
                continue
            seen.add(key)
            probes += 1
            if objective(candidate) < barrier:
                accepted.append(candidate)
                if len(accepted) >= target_count:
                    return np.asarray(accepted), probes
    if len(accepted) < 5:
        raise ValueError("could not construct a feasible optimizer population")
    return np.asarray(accepted), probes


def _string_keys(values: dict[float, float]) -> dict[str, float]:
    return {f"{prefix:.3f}": value for prefix, value in values.items()}
