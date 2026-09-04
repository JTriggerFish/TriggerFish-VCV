"""Measure and fit the production two-point frequency-dependent T60 curve."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from .decay import DecayFit, band_decay_fits, pre_onset_noise_power
from .erb import ErbFilterbank, erb_rate
from .transforms import StftConfig, stft


@dataclass(frozen=True)
class T60EnvelopeFit:
    dc_seconds: float
    nyquist_seconds: float
    log_rmse: float
    band_frequencies_hz: np.ndarray
    measured_seconds: np.ndarray
    predicted_seconds: np.ndarray
    band_fits: tuple[DecayFit, ...]


def interpolate_t60(
    frequency_hz: np.ndarray | float,
    point_frequencies_hz,
    point_seconds,
    point_active,
    sample_rate: int,
) -> np.ndarray:
    """Evaluate the C++ ERB-rate, log-T60 interpolation exactly."""
    frequencies = np.asarray(point_frequencies_hz, dtype=np.float64).copy()
    seconds = np.asarray(point_seconds, dtype=np.float64)
    active = np.asarray(point_active, dtype=bool)
    if frequencies.shape != seconds.shape or frequencies.shape != active.shape:
        raise ValueError("T60 point arrays must have matching shapes")
    if frequencies.ndim != 1 or frequencies.size < 2:
        raise ValueError("a T60 envelope needs at least two points")
    frequencies[0] = 0.0
    frequencies[-1] = 0.5 * sample_rate
    if not active[0] or not active[-1] or np.count_nonzero(active) < 2:
        raise ValueError("the DC and Nyquist T60 points must be active")
    if np.any(seconds[active] <= 0.0) or not np.isfinite(seconds[active]).all():
        raise ValueError("active T60 values must be finite and positive")
    frequencies[1:-1] = np.clip(frequencies[1:-1], 1.0, 0.5 * sample_rate - 1.0)
    seconds = np.clip(seconds, 0.01, 30.0)
    order = np.argsort(frequencies[active])
    rates = erb_rate(frequencies[active][order])
    targets = erb_rate(np.maximum(np.asarray(frequency_hz, dtype=np.float64), 0.0))
    return np.exp(np.interp(targets, rates, np.log(seconds[active][order])))


def measure_band_t60(
    samples: np.ndarray,
    sample_rate: int,
    minimum_hz: float = 60.0,
    maximum_hz: float = 20_000.0,
    band_count: int = 32,
    start_seconds: float = 0.05,
    onset_sample: int = 0,
    minimum_relative_band_db: float = -60.0,
    config: StftConfig = StftConfig(2048, 512, 4096),
) -> tuple[np.ndarray, tuple[DecayFit, ...]]:
    """Measure T60 only in perceptual bands carrying observable signal energy."""
    signal = np.asarray(samples, dtype=np.float64)
    if start_seconds < 0.0:
        raise ValueError("T60 measurement start must be non-negative")
    if not 0 <= onset_sample < signal.size:
        raise ValueError("T60 onset must lie inside the audio")
    if not np.isfinite(minimum_relative_band_db) or minimum_relative_band_db > 0.0:
        raise ValueError("minimum relative band level must be finite and non-positive")
    transform = stft(signal, sample_rate, config)
    maximum = min(maximum_hz, float(transform.frequencies_hz[-1]))
    bank = ErbFilterbank.create(
        transform.frequencies_hz, minimum_hz, maximum, band_count
    )
    energy = bank.apply_power(transform.power)
    onset_seconds = onset_sample / sample_rate
    relative_times = transform.times_seconds - onset_seconds
    half_window_seconds = 0.5 * config.window_samples / sample_rate
    pre_onset = relative_times <= -half_window_seconds
    noise = (
        pre_onset_noise_power(energy[:, pre_onset])
        if np.count_nonzero(pre_onset) >= 2
        else np.zeros(energy.shape[0])
    )
    selected = (relative_times >= start_seconds + half_window_seconds) & (
        transform.times_seconds <= signal.size / sample_rate - half_window_seconds
    )
    if np.count_nonzero(selected) < 3:
        raise ValueError("T60 measurement start leaves fewer than three frames")
    selected_energy = np.maximum(energy[:, selected] - noise[:, None], 0.0)
    band_energy = np.sum(selected_energy, axis=1)
    peak_band_energy = max(float(np.max(band_energy)), np.finfo(float).tiny)
    observable = band_energy >= peak_band_energy * 10.0 ** (
        minimum_relative_band_db / 10.0
    )
    fits = tuple(
        (
            fit
            if keep
            else replace(
                fit,
                slope_db_per_second=float("nan"),
                intercept_db=float("nan"),
                t60_seconds=float("nan"),
                r_squared=float("nan"),
                sample_count=0,
                status="below_relative_energy_floor",
            )
        )
        for fit, keep in zip(
            band_decay_fits(relative_times[selected], energy[:, selected], noise),
            observable,
        )
    )
    return bank.centers_hz, fits


def fit_two_point_t60(
    band_frequencies_hz: np.ndarray,
    band_fits: tuple[DecayFit, ...],
    sample_rate: int,
    minimum_r_squared: float = 0.95,
) -> T60EnvelopeFit:
    """Fit only DC and Nyquist log-T60 values to measured band slopes."""
    frequencies = np.asarray(band_frequencies_hz, dtype=np.float64)
    if frequencies.shape != (len(band_fits),):
        raise ValueError("T60 frequencies must match the band fits")
    measured = np.asarray([item.t60_seconds for item in band_fits])
    valid = np.asarray(
        [
            item.status == "measured"
            and np.isfinite(item.t60_seconds)
            and item.t60_seconds > 0.0
            and item.r_squared >= minimum_r_squared
            for item in band_fits
        ]
    )
    if np.count_nonzero(valid) < 2:
        raise ValueError("fewer than two reliable T60 bands were measured")

    nyquist_rate = float(erb_rate(0.5 * sample_rate))
    positions = erb_rate(frequencies[valid]) / nyquist_rate
    design = np.column_stack((1.0 - positions, positions))
    quality = np.asarray(
        [max(item.r_squared, 0.0) * max(item.sample_count, 1) for item in band_fits]
    )[valid]
    weights = np.sqrt(quality / np.max(quality))
    initial, *_ = np.linalg.lstsq(
        design * weights[:, None], np.log(measured[valid]) * weights, rcond=None
    )
    solved = least_squares(
        lambda coefficients: weights
        * (design @ coefficients - np.log(measured[valid])),
        initial,
        loss="soft_l1",
        f_scale=0.08,
    )
    coefficients = solved.x
    predicted = np.exp(design @ coefficients)
    error = np.log(predicted) - np.log(measured[valid])
    return T60EnvelopeFit(
        dc_seconds=float(np.exp(coefficients[0])),
        nyquist_seconds=float(np.exp(coefficients[1])),
        log_rmse=float(np.sqrt(np.mean(np.square(error)))),
        band_frequencies_hz=frequencies[valid],
        measured_seconds=measured[valid],
        predicted_seconds=predicted,
        band_fits=tuple(item for item, keep in zip(band_fits, valid) if keep),
    )


def recover_two_point_t60(
    samples: np.ndarray, sample_rate: int, **measurement_options
) -> T60EnvelopeFit:
    """Measure an audio decay and recover its two production T60 endpoints."""
    frequencies, fits = measure_band_t60(samples, sample_rate, **measurement_options)
    return fit_two_point_t60(frequencies, fits, sample_rate)
