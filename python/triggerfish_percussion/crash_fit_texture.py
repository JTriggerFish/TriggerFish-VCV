"""Short-time brightness and deterministic/stochastic texture comparison."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import median_filter

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer
from .crash_fit_common import fixed_rate
from .erb import ErbFilterbank
from .transforms import StftConfig, stft

_BANDS_HZ = (
    (200.0, 500.0),
    (500.0, 1000.0),
    (1000.0, 2000.0),
    (2000.0, 4000.0),
    (4000.0, 8000.0),
    (8000.0, 16000.0),
)


@dataclass(frozen=True)
class TextureQuality:
    fine_spectrum_rmse_db: float
    centroid_rmse_octaves: float
    rolloff_rmse_octaves: float
    flatness_rmse_db: float
    crest_rmse_db: float
    ridge_ratio_absolute_error: float


@dataclass(frozen=True)
class TextureDescriptor:
    fine_spectrum_db: np.ndarray
    centroid_log2_hz: np.ndarray
    rolloff_log2_hz: np.ndarray
    flatness_db: np.ndarray
    crest_db: np.ndarray
    ridge_ratio: float


def spectral_texture_quality(
    candidate: AudioBuffer,
    target: AudioBuffer,
    end_seconds: float,
    candidate_onset_sample: int | None = None,
) -> TextureQuality:
    """Compare brightness and tonal/noise texture without requiring phase identity."""
    return compare_texture_descriptors(
        texture_descriptor(candidate, end_seconds, candidate_onset_sample),
        texture_descriptor(target, end_seconds),
    )


def compare_texture_descriptors(
    candidate_features: TextureDescriptor,
    target_features: TextureDescriptor,
) -> TextureQuality:
    """Compare two reusable descriptors with no audio-side recomputation."""
    return TextureQuality(
        _salient_spectrum_rmse(
            candidate_features.fine_spectrum_db,
            target_features.fine_spectrum_db,
        ),
        _rmse(candidate_features.centroid_log2_hz - target_features.centroid_log2_hz),
        _rmse(candidate_features.rolloff_log2_hz - target_features.rolloff_log2_hz),
        _rmse(candidate_features.flatness_db - target_features.flatness_db),
        _rmse(candidate_features.crest_db - target_features.crest_db),
        abs(candidate_features.ridge_ratio - target_features.ridge_ratio),
    )


def normalized_texture_residual(
    candidate: TextureDescriptor, target: TextureDescriptor
) -> np.ndarray:
    """Return equally weighted, gate-scaled texture descriptor families."""
    groups = (
        (candidate.fine_spectrum_db - target.fine_spectrum_db) / 4.0,
        (candidate.centroid_log2_hz - target.centroid_log2_hz) / 0.35,
        (candidate.rolloff_log2_hz - target.rolloff_log2_hz) / 0.35,
        (candidate.flatness_db - target.flatness_db) / 2.5,
        (candidate.crest_db - target.crest_db) / 2.5,
        np.asarray([(candidate.ridge_ratio - target.ridge_ratio) / 0.08]),
    )
    normalized = [group / np.sqrt(max(1, group.size)) for group in groups]
    return np.concatenate(normalized) / np.sqrt(len(normalized))


def texture_descriptor(
    audio: AudioBuffer, end_seconds: float, onset_sample: int | None = None
) -> TextureDescriptor:
    """Measure one onset-aligned prefix for reuse across optimizer proposals."""
    signal, sample_rate = _onset_prefix(audio, end_seconds, onset_sample)
    config = _texture_config(signal.size)
    transform = stft(signal, sample_rate, config)
    selected = (transform.frequencies_hz >= 200.0) & (
        transform.frequencies_hz <= min(18000.0, 0.49 * sample_rate)
    )
    frequencies = transform.frequencies_hz[selected]
    power = np.maximum(transform.power[selected], np.finfo(float).tiny)
    frame_power = np.sum(power, axis=0)
    centroid = np.sum(frequencies[:, None] * power, axis=0) / frame_power
    cumulative = np.cumsum(power, axis=0) / frame_power
    rolloff = frequencies[np.argmax(cumulative >= 0.85, axis=0)]
    flatness, crest = _local_texture(power, frequencies)
    return TextureDescriptor(
        _fine_spectrum_db(transform, sample_rate),
        np.log2(np.maximum(centroid, 1.0)),
        np.log2(np.maximum(rolloff, 1.0)),
        flatness,
        crest,
        _ridge_ratio(power),
    )


def _onset_prefix(
    audio: AudioBuffer, end_seconds: float, onset_sample: int | None
) -> tuple[np.ndarray, int]:
    mono = fixed_rate(audio)
    onset = (
        detect_impact_onset(mono.samples, mono.sample_rate)
        if onset_sample is None
        else int(np.clip(onset_sample, 0, mono.samples.size))
    )
    count = max(16, round(end_seconds * mono.sample_rate))
    result = np.zeros(count, dtype=np.float64)
    available = max(0, min(count, mono.samples.size - onset))
    if available:
        result[:available] = mono.samples[onset : onset + available]
    return result, mono.sample_rate


def _texture_config(sample_count: int) -> StftConfig:
    window = 1024 if sample_count >= 1440 else (512 if sample_count >= 480 else 256)
    return StftConfig(window, window // 8, 2 * window)


def _local_texture(
    power: np.ndarray, frequencies_hz: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    flatness = []
    crest = []
    for low_hz, high_hz in _BANDS_HZ:
        selected = (frequencies_hz >= low_hz) & (frequencies_hz < high_hz)
        band = power[selected]
        if band.shape[0] < 3:
            continue
        floor = max(float(np.max(band)) * 1.0e-10, np.finfo(float).tiny)
        bounded = np.maximum(band, floor)
        arithmetic = np.mean(bounded, axis=0)
        flatness.extend(
            (10.0 * np.log10(np.exp(np.mean(np.log(bounded), axis=0)) / arithmetic))
        )
        crest.extend(10.0 * np.log10(np.max(bounded, axis=0) / arithmetic))
    return np.asarray(flatness), np.asarray(crest)


def _fine_spectrum_db(transform, sample_rate: int) -> np.ndarray:
    bank = ErbFilterbank.create(
        transform.frequencies_hz,
        200.0,
        min(16000.0, 0.49 * sample_rate),
        96,
    )
    spectrum = np.mean(bank.apply_power(transform.power), axis=1)
    level_db = 10.0 * np.log10(np.maximum(spectrum, np.finfo(float).tiny))
    return level_db - np.max(level_db)


def _salient_spectrum_rmse(candidate_db: np.ndarray, target_db: np.ndarray) -> float:
    salient = target_db >= -50.0
    difference = candidate_db[salient] - target_db[salient]
    return _rmse(difference)


def _ridge_ratio(power: np.ndarray) -> float:
    if power.shape[0] < 17 or power.shape[1] < 5:
        return 0.0
    horizontal = median_filter(power, size=(1, 7), mode="nearest")
    vertical = median_filter(power, size=(17, 1), mode="nearest")
    mask = np.square(horizontal) / (
        np.square(horizontal) + np.square(vertical) + np.finfo(float).tiny
    )
    return float(np.sum(mask * power) / np.sum(power))


def _rmse(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values)))) if values.size else 0.0
