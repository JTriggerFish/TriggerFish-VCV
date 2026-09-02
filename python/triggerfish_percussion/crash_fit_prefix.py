"""Causal, multiresolution features for cumulative crash fitting."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer
from .crash_fit_common import fixed_rate
from .crash_fit_texture import spectral_texture_quality
from .erb import ErbFilterbank
from .transforms import StftConfig, stft

CAUSAL_PREFIX_SECONDS = (0.004, 0.015, 0.100, 0.250, 1.500, 4.000)

_TIME_EDGES_SECONDS = (
    0.0,
    0.0005,
    0.001,
    0.002,
    0.004,
    0.008,
    0.015,
    0.030,
    0.060,
    0.125,
    0.250,
    0.500,
    1.000,
    1.500,
    2.500,
    4.000,
    8.000,
)

_STFT_CONFIGS = (
    StftConfig(64, 16, 128),
    StftConfig(128, 32, 256),
    StftConfig(256, 64, 512),
    StftConfig(512, 128, 1024),
    StftConfig(1024, 256, 2048),
    StftConfig(2048, 512, 4096),
)


@dataclass(frozen=True)
class CausalFeatures:
    values: np.ndarray
    weights: np.ndarray
    names: tuple[str, ...]


@dataclass(frozen=True)
class PrefixQuality:
    envelope_rmse_db: float
    spectral_rmse_db: float
    spectral_p95_absolute_db: float
    envelope_maximum_absolute_db: float = 0.0
    fine_spectrum_rmse_db: float = 0.0
    centroid_rmse_octaves: float = 0.0
    rolloff_rmse_octaves: float = 0.0
    flatness_rmse_db: float = 0.0
    crest_rmse_db: float = 0.0
    ridge_ratio_absolute_error: float = 0.0


def causal_prefix_features(
    audio: AudioBuffer, end_seconds: float, onset_sample: int | None = None
) -> CausalFeatures:
    """Describe only the onset-aligned output prefix ending at ``end_seconds``."""
    if end_seconds <= 0.0:
        raise ValueError("causal feature prefix must be positive")
    mono = fixed_rate(audio)
    onset = (
        detect_impact_onset(mono.samples, mono.sample_rate)
        if onset_sample is None
        else int(np.clip(onset_sample, 0, mono.samples.size))
    )
    signal = _prefix(mono.samples, onset, end_seconds, mono.sample_rate)
    regions = _regions(end_seconds)
    values: list[float] = []
    weights: list[float] = []
    names: list[str] = []
    _append_envelope(signal, mono.sample_rate, regions, values, weights, names)
    _append_spectra(signal, mono.sample_rate, regions, values, weights, names)
    return CausalFeatures(
        np.asarray(values, dtype=np.float64),
        np.asarray(weights, dtype=np.float64),
        tuple(names),
    )


def causal_feature_residual(
    candidate: CausalFeatures, target: CausalFeatures
) -> np.ndarray:
    if candidate.names != target.names:
        raise ValueError("causal feature layouts do not match")
    return (candidate.values - target.values) * target.weights


def causal_feature_loss(candidate: CausalFeatures, target: CausalFeatures) -> float:
    residual = causal_feature_residual(candidate, target)
    return float(np.mean(np.square(residual)))


def causal_prefix_quality(
    candidate: CausalFeatures, target: CausalFeatures
) -> PrefixQuality:
    """Return interpretable absolute errors for a fitted causal prefix."""
    if candidate.names != target.names:
        raise ValueError("causal feature layouts do not match")
    difference_db = 12.0 * (candidate.values - target.values)
    envelope = np.asarray(
        [
            value
            for value, name in zip(difference_db, target.names)
            if name.startswith("envelope/")
        ]
    )
    spectral_indices = np.asarray(
        [
            index
            for index, name in enumerate(target.names)
            if name.startswith("stft-") and target.weights[index] > 0.0
        ],
        dtype=np.intp,
    )
    spectral = difference_db[spectral_indices]
    return PrefixQuality(
        envelope_rmse_db=_rmse(envelope),
        spectral_rmse_db=_rmse(spectral),
        spectral_p95_absolute_db=(
            float(np.percentile(np.abs(spectral), 95.0)) if spectral.size else np.inf
        ),
        envelope_maximum_absolute_db=(
            float(np.max(np.abs(envelope))) if envelope.size else np.inf
        ),
    )


def causal_audio_quality(
    candidate_audio: AudioBuffer,
    target_audio: AudioBuffer,
    end_seconds: float,
    candidate_onset_sample: int | None = None,
) -> PrefixQuality:
    """Measure level/spectrum plus deterministic-versus-stochastic texture."""
    candidate = causal_prefix_features(
        candidate_audio, end_seconds, candidate_onset_sample
    )
    target = causal_prefix_features(target_audio, end_seconds)
    base = causal_prefix_quality(candidate, target)
    texture = spectral_texture_quality(
        candidate_audio,
        target_audio,
        end_seconds,
        candidate_onset_sample=candidate_onset_sample,
    )
    return PrefixQuality(
        envelope_rmse_db=base.envelope_rmse_db,
        spectral_rmse_db=base.spectral_rmse_db,
        spectral_p95_absolute_db=base.spectral_p95_absolute_db,
        envelope_maximum_absolute_db=base.envelope_maximum_absolute_db,
        fine_spectrum_rmse_db=texture.fine_spectrum_rmse_db,
        centroid_rmse_octaves=texture.centroid_rmse_octaves,
        rolloff_rmse_octaves=texture.rolloff_rmse_octaves,
        flatness_rmse_db=texture.flatness_rmse_db,
        crest_rmse_db=texture.crest_rmse_db,
        ridge_ratio_absolute_error=texture.ridge_ratio_absolute_error,
    )


def _prefix(
    samples: np.ndarray, onset: int, end_seconds: float, sample_rate: int
) -> np.ndarray:
    count = max(8, round(end_seconds * sample_rate))
    result = np.zeros(count, dtype=np.float64)
    available = max(0, min(count, samples.size - onset))
    if available:
        result[:available] = samples[onset : onset + available]
    return result


def _regions(end_seconds: float) -> tuple[tuple[float, float], ...]:
    clipped = [value for value in _TIME_EDGES_SECONDS if value < end_seconds]
    clipped.append(end_seconds)
    return tuple(
        (start, end) for start, end in zip(clipped, clipped[1:]) if end > start
    )


def _append_envelope(
    signal: np.ndarray,
    sample_rate: int,
    regions: tuple[tuple[float, float], ...],
    values: list[float],
    weights: list[float],
    names: list[str],
) -> None:
    normalization = 1.0 / np.sqrt(max(1, 2 * len(regions)))
    floor = 1.0e-12
    for start, end in regions:
        first = min(signal.size, round(start * sample_rate))
        last = min(signal.size, max(first + 1, round(end * sample_rate)))
        frame = signal[first:last]
        rms = np.sqrt(np.mean(np.square(frame)))
        peak = np.max(np.abs(frame))
        crest = peak / max(rms, floor)
        label = f"{1000.0 * start:.3f}-{1000.0 * end:.3f}ms"
        values.extend(
            (
                20.0 * np.log10(max(rms, floor)) / 12.0,
                np.log(max(crest, floor)),
            )
        )
        weights.extend((normalization, 0.35 * normalization))
        names.extend((f"envelope/{label}", f"crest/{label}"))


def _append_spectra(
    signal: np.ndarray,
    sample_rate: int,
    regions: tuple[tuple[float, float], ...],
    values: list[float],
    weights: list[float],
    names: list[str],
) -> None:
    configs = _configs_for_prefix(signal.size)
    feature_count = max(1, len(configs) * len(regions) * 32)
    normalization = 1.0 / np.sqrt(feature_count)
    for config in configs:
        transform = stft(signal, sample_rate, config)
        bank = ErbFilterbank.create(
            transform.frequencies_hz,
            80.0,
            min(18000.0, 0.49 * sample_rate),
            32,
        )
        band_power = bank.apply_power(transform.power)
        event_peak = max(float(np.max(band_power)), 1.0e-24)
        for start, end in regions:
            selected = (transform.times_seconds >= start) & (
                transform.times_seconds < end
            )
            if np.any(selected):
                spectrum = np.mean(band_power[:, selected], axis=1)
            else:
                frame = int(np.argmin(np.abs(transform.times_seconds - start)))
                spectrum = band_power[:, frame]
            level = 10.0 * np.log10(np.maximum(spectrum, 1.0e-24)) / 12.0
            relative_db = 10.0 * np.log10(np.maximum(spectrum, 1.0e-24) / event_peak)
            salience = np.clip((relative_db + 90.0) / 40.0, 0.0, 1.0)
            label = f"{1000.0 * start:.3f}-{1000.0 * end:.3f}ms"
            values.extend(level.tolist())
            weights.extend((normalization * salience).tolist())
            names.extend(
                f"stft-{config.window_samples}/{label}/{frequency:.0f}Hz"
                for frequency in bank.centers_hz
            )


def _configs_for_prefix(sample_count: int) -> tuple[StftConfig, ...]:
    maximum_window = max(128, min(2048, 2 * sample_count))
    selected = tuple(
        config for config in _STFT_CONFIGS if config.window_samples <= maximum_window
    )
    return selected or (_STFT_CONFIGS[0],)


def _rmse(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values)))) if values.size else np.inf
