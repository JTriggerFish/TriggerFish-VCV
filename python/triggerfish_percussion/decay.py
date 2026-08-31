"""Noise-aware energy-decay curves and robust decay fits."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.stats import linregress, theilslopes


@dataclass(frozen=True)
class DecayFit:
    slope_db_per_second: float
    intercept_db: float
    t60_seconds: float
    r_squared: float
    sample_count: int
    status: str = "measured"
    fit_end_seconds: float = float("nan")
    noise_limit_seconds: float = float("nan")
    tail_correction_fraction: float = float("nan")


def energy_decay_db(energy: np.ndarray, noise_power: float = 0.0) -> np.ndarray:
    values = np.asarray(energy, dtype=np.float64)
    if (
        values.ndim != 1
        or values.size < 2
        or np.any(values < 0)
        or not np.isfinite(values).all()
    ):
        raise ValueError("decay input must be finite non-negative energy")
    noise = float(noise_power)
    if not np.isfinite(noise) or noise < 0:
        raise ValueError("decay noise power must be finite and non-negative")
    corrected = np.maximum(values - noise, 0.0)
    integrated = np.cumsum(corrected[::-1])[::-1]
    peak = max(float(integrated[0]), np.finfo(float).tiny)
    floor = max(peak * 1.0e-12, np.finfo(float).tiny)
    return 10.0 * np.log10(np.maximum(integrated, floor) / peak)


def fit_decay(
    times_seconds: np.ndarray,
    decay_db: np.ndarray,
    upper_db: float = -5.0,
    lower_db: float = -25.0,
    robust: bool = True,
) -> DecayFit:
    times = np.asarray(times_seconds, dtype=np.float64)
    decay = np.asarray(decay_db, dtype=np.float64)
    if (
        times.shape != decay.shape
        or times.ndim != 1
        or not np.isfinite(times).all()
        or not np.isfinite(decay).all()
        or not np.all(np.diff(times) > 0)
    ):
        raise ValueError("decay times must increase and match finite levels")
    if not lower_db < upper_db:
        raise ValueError("decay fit requires lower_db below upper_db")
    selected = (decay <= upper_db) & (decay >= lower_db)
    if np.count_nonzero(selected) < 3:
        return _invalid_fit("insufficient_dynamic_range", np.count_nonzero(selected))
    x = times[selected]
    y = decay[selected]
    if robust:
        slope, intercept, _, _ = theilslopes(y, x)
    else:
        regression = linregress(x, y)
        slope, intercept = regression.slope, regression.intercept
    prediction = slope * x + intercept
    residual = float(np.sum(np.square(y - prediction)))
    total = float(np.sum(np.square(y - np.mean(y))))
    r_squared = 1.0 - residual / total if total > 0 else 1.0
    t60 = -60.0 / slope if slope < 0 else float("inf")
    if not np.isfinite(t60):
        return _invalid_fit("non_decaying", x.size)
    return DecayFit(
        float(slope),
        float(intercept),
        float(t60),
        r_squared,
        x.size,
        fit_end_seconds=float(x[-1]),
    )


def pre_onset_noise_power(band_energy: np.ndarray) -> np.ndarray:
    """Return one robust noise-power estimate per band from pre-onset frames."""
    energy = np.asarray(band_energy, dtype=np.float64)
    if (
        energy.ndim != 2
        or energy.shape[1] < 2
        or np.any(energy < 0)
        or not np.isfinite(energy).all()
    ):
        raise ValueError("pre-onset noise requires finite band-by-time energy")
    return np.median(energy, axis=1)


def band_decay_fits(
    times_seconds: np.ndarray,
    band_energy: np.ndarray,
    noise_power: np.ndarray | None = None,
    upper_db: float = -5.0,
    lower_db: float = -25.0,
) -> tuple[DecayFit, ...]:
    times = np.asarray(times_seconds, dtype=np.float64)
    energy = np.asarray(band_energy, dtype=np.float64)
    if energy.ndim != 2 or energy.shape[1] != times.size:
        raise ValueError("band decay expects band-by-time energy")
    if np.any(energy < 0) or not np.isfinite(energy).all():
        raise ValueError("band decay energy must be finite and non-negative")
    noise = (
        np.zeros(energy.shape[0], dtype=np.float64)
        if noise_power is None
        else np.asarray(noise_power, dtype=np.float64)
    )
    if noise.shape != (energy.shape[0],):
        raise ValueError("band noise power must provide one value per band")
    if np.any(noise < 0) or not np.isfinite(noise).all():
        raise ValueError("band noise power must be finite and non-negative")
    return tuple(
        _band_decay_fit(times, band, floor, upper_db, lower_db)
        for band, floor in zip(energy, noise)
    )


def _band_decay_fit(
    times: np.ndarray,
    energy: np.ndarray,
    noise_power: float,
    upper_db: float,
    lower_db: float,
) -> DecayFit:
    if times.size < 3 or not np.all(np.diff(times) > 0):
        raise ValueError("band decay times must increase")
    corrected = np.maximum(energy - noise_power, 0.0)
    integration_end = _integration_end(energy, noise_power)
    intersection = float(times[integration_end]) if noise_power > 0 else float("nan")
    if integration_end < 2:
        return replace(
            _invalid_fit("noise_limited", 0),
            noise_limit_seconds=float(times[integration_end]),
        )
    tail_energy = _extrapolated_tail_energy(
        times[: integration_end + 1],
        corrected[: integration_end + 1],
        9.0 * noise_power,
    )
    integrated = np.cumsum(corrected[: integration_end + 1][::-1])[::-1]
    integrated += tail_energy
    if integrated[0] <= np.finfo(float).tiny:
        return replace(
            _invalid_fit("insufficient_energy", 0),
            noise_limit_seconds=intersection,
        )
    peak = max(float(integrated[0]), np.finfo(float).tiny)
    decay = 10.0 * np.log10(np.maximum(integrated, peak * 1.0e-12) / peak)
    reached_interval = float(np.min(decay)) <= lower_db
    if not reached_interval:
        status = "noise_limited" if noise_power > 0 else "recording_too_short"
        return replace(
            _invalid_fit(status, np.count_nonzero(decay <= upper_db)),
            noise_limit_seconds=intersection,
        )
    fitted = fit_decay(times[: integration_end + 1], decay, upper_db, lower_db)
    return replace(
        fitted,
        noise_limit_seconds=intersection,
        tail_correction_fraction=float(tail_energy / integrated[0]),
    )


def _integration_end(energy: np.ndarray, noise_power: float) -> int:
    if noise_power <= 0:
        return energy.size - 1
    above = np.flatnonzero(energy >= 10.0 * noise_power)
    return int(above[-1]) if above.size else 0


def _extrapolated_tail_energy(
    times: np.ndarray, energy: np.ndarray, minimum_energy: float
) -> float:
    positive = np.flatnonzero(energy > minimum_energy)
    if positive.size < 8:
        return 0.0
    count = min(positive.size, max(32, positive.size // 4))
    selected = positive[-count:]
    levels = 10.0 * np.log10(energy[selected])
    slope, intercept, _, _ = theilslopes(levels, times[selected])
    step = float(np.median(np.diff(times)))
    if not np.isfinite(slope) or slope >= 0 or step <= 0:
        return 0.0
    ratio = 10.0 ** (float(slope) * step / 10.0)
    if not 0.0 < ratio < 1.0:
        return 0.0
    predicted = 10.0 ** ((float(intercept) + float(slope) * times[-1]) / 10.0)
    tail = predicted * ratio / (1.0 - ratio)
    return float(tail) if np.isfinite(tail) else 0.0


def _invalid_fit(status: str, sample_count: int) -> DecayFit:
    return DecayFit(
        float("nan"),
        float("nan"),
        float("nan"),
        float("nan"),
        int(sample_count),
        status,
    )
