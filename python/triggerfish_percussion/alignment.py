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
    minimum_peak_energy_ratio: float = 1.0e-5,
) -> int:
    """Locate the first perceptually material impact energy."""
    signal = _mono_signal(samples)
    if not 0.0 < minimum_peak_energy_ratio < 1.0:
        raise ValueError("impact peak-energy ratio must lie between zero and one")
    width = max(1, int(round(smoothing_seconds * sample_rate)))
    window = np.ones(width) / width
    full_energy = np.convolve(np.square(signal), window, mode="same")
    full_peak = float(np.max(full_energy))
    if full_peak <= np.finfo(float).tiny:
        return 0
    # A novelty-only detector locked onto inaudible edit noise in several sample
    # layers, putting the fitted strike 6--9 ms before the actual contact. The
    # -50 dB power criterion is low enough to retain a quiet cymbal contact while
    # rejecting such isolated pre-roll artifacts.
    active = full_energy >= minimum_peak_energy_ratio * full_peak
    candidates = np.flatnonzero(active)
    coarse = int(candidates[0]) if candidates.size else int(np.argmax(full_energy))

    # Refine against high-frequency contact energy, but normalize only inside a
    # short interval after the coarse onset. A much brighter later bloom then
    # cannot move the detector, while low-frequency pre-contact motion does not
    # become the cymbal strike anchor.
    cutoff = min(highpass_hz, 0.4 * sample_rate)
    filtered = sosfilt(
        butter(2, cutoff, btype="highpass", fs=sample_rate, output="sos"), signal
    )
    high_energy = np.convolve(np.square(filtered), window, mode="same")
    local_first = coarse
    local_last = min(signal.size, coarse + round(0.010 * sample_rate))
    local_peak = float(np.max(high_energy[local_first:local_last]))
    local = np.flatnonzero(
        high_energy[local_first:local_last] >= minimum_peak_energy_ratio * local_peak
    )
    onset = local_first + int(local[0]) if local.size else coarse
    # ``mode='same'`` centers the smoothing window, so a nonzero detection is
    # early by roughly half a window. Do not compensate an event that genuinely
    # starts at the first available sample: there is no discarded pre-roll.
    return min(signal.size - 1, onset + width // 2) if onset > 0 else 0


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
