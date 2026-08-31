"""Onset detection and sample-preserving pair alignment."""

from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, correlate, correlation_lags, sosfilt


@dataclass(frozen=True)
class Alignment:
    reference_onset: int
    candidate_onset: int
    candidate_lag_samples: int


def detect_impact_onset(
    samples: np.ndarray,
    sample_rate: float,
    highpass_hz: float = 1200.0,
    smoothing_seconds: float = 0.00075,
) -> int:
    signal = _mono_signal(samples)
    cutoff = min(highpass_hz, 0.4 * sample_rate)
    filtered = sosfilt(
        butter(2, cutoff, btype="highpass", fs=sample_rate, output="sos"), signal
    )
    width = max(1, int(round(smoothing_seconds * sample_rate)))
    energy = np.convolve(np.square(filtered), np.ones(width) / width, mode="same")
    novelty = np.maximum(
        np.diff(np.sqrt(energy + np.finfo(float).tiny), prepend=0.0), 0.0
    )
    peak = float(np.max(novelty))
    if peak <= np.finfo(float).tiny:
        return 0
    baseline_count = max(width, min(signal.size // 10, int(0.05 * sample_rate)))
    baseline = novelty[:baseline_count]
    median = float(np.median(baseline))
    deviation = float(np.median(np.abs(baseline - median)))
    # Cymbals can bloom much louder than their first contact. A high fraction
    # of the global peak therefore skips the physical onset; the pre-onset
    # robust floor is primary and the peak term only rejects numerical dust.
    threshold = max(median + 8.0 * deviation, 1.0e-6 * peak)
    candidates = np.flatnonzero(novelty >= threshold)
    return int(candidates[0]) if candidates.size else int(np.argmax(novelty))


def measure_alignment(
    reference: np.ndarray,
    candidate: np.ndarray,
    sample_rate: float,
    correlation_seconds: float = 0.02,
    maximum_lag_seconds: float = 0.005,
    refinement: str = "onset",
) -> Alignment:
    reference_signal = _mono_signal(reference)
    candidate_signal = _mono_signal(candidate)
    reference_onset = detect_impact_onset(reference_signal, sample_rate)
    candidate_onset = detect_impact_onset(candidate_signal, sample_rate)
    onset_lag = reference_onset - candidate_onset
    if refinement == "onset":
        return Alignment(reference_onset, candidate_onset, onset_lag)
    if refinement != "waveform_copy":
        raise ValueError(f"unknown alignment refinement: {refinement}")
    radius = max(8, int(round(correlation_seconds * sample_rate)))
    maximum_lag = max(1, int(round(maximum_lag_seconds * sample_rate)))
    reference_window = _window(reference_signal, reference_onset, radius)
    candidate_window = _window(candidate_signal, candidate_onset, radius)
    correlation = correlate(
        reference_window, candidate_window, mode="full", method="fft"
    )
    lags = correlation_lags(reference_window.size, candidate_window.size, mode="full")
    selected = np.abs(lags) <= maximum_lag
    local_lag = int(lags[selected][np.argmax(correlation[selected])])
    candidate_lag = onset_lag - local_lag
    return Alignment(reference_onset, candidate_onset, candidate_lag)


def shift_with_zeros(samples: np.ndarray, lag_samples: int) -> np.ndarray:
    signal = np.asarray(samples, dtype=np.float64)
    shifted = np.zeros_like(signal)
    if abs(lag_samples) >= signal.shape[0]:
        return shifted
    if lag_samples >= 0:
        shifted[lag_samples:] = signal[: signal.shape[0] - lag_samples]
    else:
        shifted[:lag_samples] = signal[-lag_samples:]
    return shifted


def _mono_signal(samples: np.ndarray) -> np.ndarray:
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or signal.size < 2 or not np.isfinite(signal).all():
        raise ValueError("alignment requires finite non-empty mono audio")
    return signal


def _window(signal: np.ndarray, center: int, radius: int) -> np.ndarray:
    output = np.zeros(2 * radius, dtype=np.float64)
    source_start = max(0, center - radius // 4)
    source_end = min(signal.size, source_start + output.size)
    output[: source_end - source_start] = signal[source_start:source_end]
    return output
