"""Central finite-difference conditioning audit for crash fitting controls."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .crash_fit_common import CrashFitCell, render_cell
from .crash_fit_parameters import (
    ATTACK_PARAMETERS,
    FitParameter,
    fit_parameter_value,
    replace_fit_parameters,
)
from .crash_fit_prefix import (
    causal_feature_residual,
    causal_prefix_quality,
    causal_prefix_features,
)
from .crash_fit_texture import (
    compare_texture_descriptors,
    normalized_texture_residual,
    texture_descriptor,
)
from .crash_model import CrashFit

AUDIT_PREFIX_SECONDS = (0.004, 0.015, 0.030, 0.100)


def parameter_conditioning_diagnostics(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    parameters: tuple[FitParameter, ...] = ATTACK_PARAMETERS,
    prefixes: tuple[float, ...] = AUDIT_PREFIX_SECONDS,
) -> dict[str, object]:
    """Measure range-normalized feature Jacobians and correlated directions."""
    if not cells or not parameters or not prefixes:
        raise ValueError("conditioning audit inputs cannot be empty")
    targets = _targets(cells, prefixes)
    derivatives = []
    results = []
    for parameter in parameters:
        low, high, range_fraction = _central_points(fit, parameter, 1.0)
        low_vectors = _residual_vectors(cells, low, prefixes, targets)
        high_vectors = _residual_vectors(cells, high, prefixes, targets)
        prefix_norms = {}
        parameter_derivatives = []
        loss_slopes = {}
        for prefix in prefixes:
            derivative = (high_vectors[prefix] - low_vectors[prefix]) / range_fraction
            parameter_derivatives.append(derivative)
            prefix_norms[f"{prefix:.3f}"] = _rms(derivative)
            loss_slopes[f"{prefix:.3f}"] = (
                _mean_square(high_vectors[prefix]) - _mean_square(low_vectors[prefix])
            ) / range_fraction
        combined = np.concatenate(parameter_derivatives)
        fine_low, fine_high, fine_fraction = _central_points(fit, parameter, 0.5)
        fine_low_vectors = _residual_vectors(cells, fine_low, prefixes, targets)
        fine_high_vectors = _residual_vectors(cells, fine_high, prefixes, targets)
        fine_combined = np.concatenate(
            [
                (fine_high_vectors[prefix] - fine_low_vectors[prefix]) / fine_fraction
                for prefix in prefixes
            ]
        )
        derivatives.append(combined)
        results.append(
            {
                "parameter": parameter.name,
                "low": fit_parameter_value(low, parameter.name),
                "high": fit_parameter_value(high, parameter.name),
                "range_fraction": range_fraction,
                "jacobian_rms": _rms(combined),
                "half_step_jacobian_rms": _rms(fine_combined),
                "step_stability_cosine": _cosine(combined, fine_combined),
                "step_norm_delta_db": float(
                    20.0
                    * np.log10(
                        max(_rms(fine_combined), 1.0e-30) / max(_rms(combined), 1.0e-30)
                    )
                ),
                "prefix_jacobian_rms": prefix_norms,
                "prefix_loss_slopes": loss_slopes,
            }
        )
    maximum = max(item["jacobian_rms"] for item in results)
    for item in results:
        item["relative_leverage_db"] = float(
            20.0 * np.log10(max(item["jacobian_rms"], 1.0e-30) / maximum)
        )
    return {
        "cell_count": len(cells),
        "prefix_seconds": list(prefixes),
        "finite_difference": "central secant, normalized by declared range",
        "parameters": results,
        "strongest_correlations": _strongest_correlations(
            parameters, derivatives, results
        ),
    }


def morris_screening_diagnostics(
    cells: tuple[CrashFitCell, ...],
    fit: CrashFit,
    parameters: tuple[FitParameter, ...] = ATTACK_PARAMETERS,
    prefixes: tuple[float, ...] = AUDIT_PREFIX_SECONDS,
    trajectories: int = 8,
    levels: int = 6,
    seed: int = 20260901,
    workers: int = 8,
) -> dict[str, object]:
    """Globally screen nonlinear controls with optimized Morris trajectories."""
    try:
        from SALib.analyze import morris as morris_analyze
        from SALib.sample import morris as morris_sample
    except ImportError as error:
        raise RuntimeError(
            "install the 'analysis' extra to run Morris screening"
        ) from error
    if trajectories < 2 or levels < 4 or levels % 2 or workers < 1:
        raise ValueError("Morris needs >=2 trajectories, even levels >=4, and workers")
    problem = {
        "num_vars": len(parameters),
        "names": [parameter.name for parameter in parameters],
        "bounds": [[parameter.lower, parameter.upper] for parameter in parameters],
    }
    pool_size = max(4 * trajectories, trajectories + 1)
    samples = morris_sample.sample(
        problem,
        pool_size,
        num_levels=levels,
        optimal_trajectories=trajectories,
        local_optimization=True,
        seed=seed,
    )
    targets = _targets(cells, prefixes)

    def evaluate(values):
        candidate = replace_fit_parameters(
            fit,
            {
                parameter.name: float(value)
                for parameter, value in zip(parameters, values)
            },
        )
        scores = _screening_scores(cells, candidate, prefixes, targets)
        return [
            scores[prefix][family] for prefix in prefixes for family in scores[prefix]
        ]

    with ThreadPoolExecutor(max_workers=workers) as executor:
        outputs = np.asarray(list(executor.map(evaluate, samples)))
    results = {}
    output_index = 0
    family_names = _screening_family_names(prefixes)
    for prefix, families in zip(prefixes, family_names):
        prefix_results = {}
        for family in families:
            analysis = morris_analyze.analyze(
                problem,
                samples,
                outputs[:, output_index],
                num_resamples=200,
                conf_level=0.95,
                scaled=True,
                num_levels=levels,
                seed=seed + output_index,
            )
            prefix_results[family] = [
                {
                    "parameter": name,
                    "mu": float(analysis["mu"][parameter_index]),
                    "mu_star": float(analysis["mu_star"][parameter_index]),
                    "sigma": float(analysis["sigma"][parameter_index]),
                    "mu_star_confidence_95": float(
                        analysis["mu_star_conf"][parameter_index]
                    ),
                }
                for parameter_index, name in enumerate(problem["names"])
            ]
            output_index += 1
        results[f"{prefix:.3f}"] = prefix_results
    return {
        "method": "SALib optimized Morris elementary effects",
        "candidate_evaluations": int(samples.shape[0]),
        "trajectories": trajectories,
        "levels": levels,
        "workers": workers,
        "prefix_seconds": list(prefixes),
        "results": results,
    }


def _screening_scores(cells, fit, prefixes, targets):
    """Return gate-normalized loss families from one render per cell."""
    scores = {
        prefix: {family: [] for family in _families_for_prefix(prefix)}
        for prefix in prefixes
    }
    duration = max(prefixes) + 0.05
    for index, cell in enumerate(cells):
        rendered = render_cell(cell, fit, duration)
        for prefix in prefixes:
            target_causal, target_texture = targets[(index, prefix)]
            features = causal_prefix_features(rendered, prefix, onset_sample=0)
            residual = causal_feature_residual(features, target_causal)
            quality = causal_prefix_quality(features, target_causal)
            scores[prefix]["causal"].append(_mean_square(residual))
            scores[prefix]["envelope"].append(np.square(quality.envelope_rmse_db / 6.0))
            scores[prefix]["coarse_spectrum"].append(
                np.mean(
                    np.square(
                        (
                            quality.spectral_rmse_db / 10.0,
                            quality.spectral_p95_absolute_db / 18.0,
                        )
                    )
                )
            )
            if target_texture is None:
                continue
            texture = compare_texture_descriptors(
                texture_descriptor(rendered, prefix, onset_sample=0),
                target_texture,
            )
            scores[prefix]["brightness"].append(
                np.mean(
                    np.square(
                        (
                            texture.fine_spectrum_rmse_db / 4.0,
                            texture.centroid_rmse_octaves / 0.35,
                            texture.rolloff_rmse_octaves / 0.35,
                        )
                    )
                )
            )
            scores[prefix]["texture"].append(
                np.mean(
                    np.square(
                        (
                            texture.flatness_rmse_db / 2.5,
                            texture.crest_rmse_db / 2.5,
                            texture.ridge_ratio_absolute_error / 0.08,
                        )
                    )
                )
            )
    return {
        prefix: {family: float(np.mean(values)) for family, values in families.items()}
        for prefix, families in scores.items()
    }


def _screening_family_names(prefixes):
    return tuple(_families_for_prefix(prefix) for prefix in prefixes)


def _families_for_prefix(prefix):
    base = ("causal", "envelope", "coarse_spectrum")
    return base + (("brightness", "texture") if prefix >= 0.015 else ())


def _targets(cells, prefixes):
    return {
        (index, prefix): (
            causal_prefix_features(cell.reference, prefix),
            texture_descriptor(cell.reference, prefix) if prefix >= 0.015 else None,
        )
        for index, cell in enumerate(cells)
        for prefix in prefixes
    }


def _residual_vectors(cells, fit, prefixes, targets):
    by_prefix = {prefix: [] for prefix in prefixes}
    duration = max(prefixes) + 0.05
    for index, cell in enumerate(cells):
        rendered = render_cell(cell, fit, duration)
        for prefix in prefixes:
            target_causal, target_texture = targets[(index, prefix)]
            causal = causal_feature_residual(
                causal_prefix_features(rendered, prefix, onset_sample=0),
                target_causal,
            )
            parts = [causal]
            if target_texture is not None:
                parts.append(
                    normalized_texture_residual(
                        texture_descriptor(rendered, prefix, onset_sample=0),
                        target_texture,
                    )
                )
            by_prefix[prefix].append(np.concatenate(parts))
    return {prefix: np.concatenate(values) for prefix, values in by_prefix.items()}


def _central_points(
    fit: CrashFit, parameter: FitParameter, step_scale: float
) -> tuple[CrashFit, CrashFit, float]:
    value = fit_parameter_value(fit, parameter.name)
    span = parameter.upper - parameter.lower
    step = max(step_scale * parameter.relative_step * span, 1.0e-4 * span)
    low_value = max(parameter.lower, value - step)
    high_value = min(parameter.upper, value + step)
    if high_value <= low_value:
        raise ValueError(f"parameter {parameter.name} has no finite-difference span")
    return (
        replace_fit_parameters(fit, {parameter.name: low_value}),
        replace_fit_parameters(fit, {parameter.name: high_value}),
        (high_value - low_value) / span,
    )


def _strongest_correlations(parameters, derivatives, results):
    pairs = []
    for first in range(len(parameters)):
        if results[first]["relative_leverage_db"] < -40.0:
            continue
        for second in range(first + 1, len(parameters)):
            if results[second]["relative_leverage_db"] < -40.0:
                continue
            cosine = _cosine(derivatives[first], derivatives[second])
            pairs.append(
                {
                    "first": parameters[first].name,
                    "second": parameters[second].name,
                    "cosine": cosine,
                }
            )
    return sorted(pairs, key=lambda item: abs(item["cosine"]), reverse=True)[:30]


def _cosine(first: np.ndarray, second: np.ndarray) -> float:
    denominator = np.linalg.norm(first) * np.linalg.norm(second)
    return float(np.dot(first, second) / denominator) if denominator > 0.0 else 0.0


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def _mean_square(values: np.ndarray) -> float:
    return float(np.mean(np.square(values)))
