"""Loss construction and hard gates for crash object-profile refinement."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .crash_fit_common import CrashFitCell, render_cell
from .crash_fit_parameters import (
    BODY_DECAY_SLOTS,
    fit_parameter_value,
    replace_fit_parameters,
)
from .crash_fit_prefix import (
    PrefixQuality,
    causal_audio_quality,
    causal_feature_loss,
    causal_feature_residual,
    causal_prefix_features,
)
from .crash_fit_texture import (
    normalized_texture_residual,
    texture_descriptor,
)
from .crash_model import CrashFit

PREFIXES = (0.004, 0.015, 0.100)
PROTECTED_LOSS_TOLERANCE = 0.01


@dataclass(frozen=True)
class BoundedParameter:
    name: str
    lower: float
    upper: float


TEMPORAL_PARAMETERS = (
    *(
        BoundedParameter(f"body_decay_seconds[{index}]", 0.02, 20.0)
        for index in BODY_DECAY_SLOTS
    ),
    *(BoundedParameter(f"turbulence_gain[{index}]", 0.0, 1.0) for index in range(3)),
    BoundedParameter("turbulence_persistence", 0.25, 4.0),
    BoundedParameter("colour_frequency_hz", 1000.0, 14000.0),
    BoundedParameter("colour_gain_db", -12.0, 12.0),
    BoundedParameter("high_cut_hz", 6000.0, 22000.0),
)


@dataclass(frozen=True)
class SpectralTargets:
    causal: dict[float, object]
    texture: dict[float, object]


def temporal_spectral_parameter_names() -> tuple[str, ...]:
    """Expose the constrained second-pass controls for audit tests."""
    return tuple(parameter.name for parameter in TEMPORAL_PARAMETERS)


def temporal_bounds(initial: CrashFit):
    names = temporal_spectral_parameter_names()
    lower = np.asarray([item.lower for item in TEMPORAL_PARAMETERS])
    upper = np.asarray([item.upper for item in TEMPORAL_PARAMETERS])
    start = np.asarray([fit_parameter_value(initial, name) for name in names])
    if np.any(start < lower) or np.any(start > upper):
        raise ValueError("temporal spectrum seed lies outside its declared bounds")
    return names, lower, upper, start


def replace_temporal_parameters(initial: CrashFit, names, values):
    return replace_fit_parameters(initial, dict(zip(names, values.tolist())))


def make_spectral_targets(cell: CrashFitCell) -> SpectralTargets:
    return SpectralTargets(
        {prefix: causal_prefix_features(cell.reference, prefix) for prefix in PREFIXES},
        {
            prefix: texture_descriptor(cell.reference, prefix)
            for prefix in PREFIXES
            if prefix >= 0.015
        },
    )


def profile_residual(cell, fit, targets, start):
    audio = render_cell(cell, fit, 0.15)
    parts = _causal_residuals(
        audio, targets, ((0.004, 0.25), (0.015, 0.6), (0.1, 1.5)), 1.0
    )
    parts += _texture_residuals(audio, targets, ((0.015, 0.5), (0.1, 2.0)), 1.25, 0.0)
    envelope = np.asarray(fit.dense_gain_envelope_db)
    parts.append(0.18 * np.diff(envelope, n=2) / (6.0 * np.sqrt(envelope.size - 2)))
    parts.append(0.04 * (envelope[1:] - start) / (12.0 * np.sqrt(start.size)))
    return np.concatenate(parts)


def temporal_residual(cell, fit, targets, values, start, lower, upper):
    audio = render_cell(cell, fit, 0.15)
    parts = _causal_residuals(
        audio, targets, ((0.004, 0.25), (0.015, 0.6), (0.1, 2.2)), 1.5
    )
    parts += _texture_residuals(audio, targets, ((0.015, 0.6), (0.1, 2.0)), 1.0, 1.5)
    scale = np.maximum(upper - lower, np.finfo(float).tiny)
    parts.append(0.01 * (values - start) / (scale * np.sqrt(values.size)))
    return np.concatenate(parts)


def prefix_qualities(cell, fit):
    audio = render_cell(cell, fit, 0.15)
    return {
        prefix: causal_audio_quality(
            audio, cell.reference, prefix, candidate_onset_sample=0
        )
        for prefix in PREFIXES
    }


def quality_passes(quality: PrefixQuality, prefix: float) -> bool:
    base = (
        quality.envelope_rmse_db <= 6.0
        and quality.envelope_maximum_absolute_db <= 6.0
        and quality.spectral_rmse_db <= 10.0
        and quality.spectral_p95_absolute_db <= 18.0
    )
    return base and (
        prefix < 0.015
        or (
            quality.fine_spectrum_rmse_db <= 4.0
            and quality.centroid_rmse_octaves <= 0.35
            and quality.rolloff_rmse_octaves <= 0.35
            and quality.flatness_rmse_db <= 2.5
            and quality.crest_rmse_db <= 2.5
            and quality.ridge_ratio_absolute_error <= 0.08
        )
    )


def prefix_losses(cell, fit, targets):
    audio = render_cell(cell, fit, 0.15)
    return {
        prefix: causal_feature_loss(
            causal_prefix_features(audio, prefix, onset_sample=0),
            targets.causal[prefix],
        )
        for prefix in PREFIXES
    }


def _causal_residuals(audio, targets, weighted_prefixes, envelope_weight):
    parts = []
    for prefix, weight in weighted_prefixes:
        candidate = causal_prefix_features(audio, prefix, onset_sample=0)
        target = targets.causal[prefix]
        parts.append(weight * causal_feature_residual(candidate, target))
        envelope = np.asarray(
            [
                2.0 * (value - expected)
                for value, expected, name in zip(
                    candidate.values, target.values, target.names
                )
                if name.startswith("envelope/")
            ]
        )
        parts.append(envelope_weight * weight * envelope / np.sqrt(envelope.size))
    return parts


def _texture_residuals(audio, targets, weighted_prefixes, fine_weight, centroid_weight):
    parts = []
    for prefix, weight in weighted_prefixes:
        candidate = texture_descriptor(audio, prefix, onset_sample=0)
        target = targets.texture[prefix]
        parts.append(weight * normalized_texture_residual(candidate, target))
        salient = target.fine_spectrum_db >= -50.0
        fine = candidate.fine_spectrum_db[salient] - target.fine_spectrum_db[salient]
        parts.append(fine_weight * weight * fine / (4.0 * np.sqrt(fine.size)))
        if centroid_weight:
            centroid = candidate.centroid_log2_hz - target.centroid_log2_hz
            parts.append(
                centroid_weight * weight * centroid / (0.35 * np.sqrt(centroid.size))
            )
    return parts
