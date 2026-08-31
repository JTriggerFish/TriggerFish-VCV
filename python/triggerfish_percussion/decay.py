"""Noise-aware energy-decay curves and robust decay fits."""

from dataclasses import dataclass

import numpy as np
from scipy.stats import linregress, theilslopes


@dataclass(frozen=True)
class DecayFit:
    slope_db_per_second: float
    intercept_db: float
    t60_seconds: float
    r_squared: float
    sample_count: int


def energy_decay_db(energy: np.ndarray, noise_power: float = 0.0) -> np.ndarray:
    values = np.asarray(energy, dtype=np.float64)
    if (
        values.ndim != 1
        or values.size < 2
        or np.any(values < 0)
        or not np.isfinite(values).all()
    ):
        raise ValueError("decay input must be finite non-negative energy")
    noise = max(0.0, float(noise_power))
    corrected = np.maximum(values - noise, 0.0)
    integrated = np.cumsum(corrected[::-1])[::-1]
    peak = max(float(integrated[0]), np.finfo(float).tiny)
    floor = peak * 1.0e-12
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
    ):
        raise ValueError("decay time and level vectors must match and be finite")
    selected = (decay <= upper_db) & (decay >= lower_db)
    if np.count_nonzero(selected) < 3:
        return DecayFit(
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
            int(np.count_nonzero(selected)),
        )
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
    return DecayFit(float(slope), float(intercept), float(t60), r_squared, x.size)


def estimate_noise_power(energy: np.ndarray, fraction: float = 0.1) -> float:
    values = np.asarray(energy, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or np.any(values < 0):
        raise ValueError("noise estimate requires non-negative energy")
    count = max(1, int(round(np.clip(fraction, 0.01, 0.5) * values.size)))
    return float(np.median(values[-count:]))


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
    noise = (
        np.array([estimate_noise_power(band) for band in energy])
        if noise_power is None
        else np.asarray(noise_power, dtype=np.float64)
    )
    if noise.shape != (energy.shape[0],):
        raise ValueError("band noise power must provide one value per band")
    return tuple(
        fit_decay(times, energy_decay_db(band, floor), upper_db, lower_db)
        for band, floor in zip(energy, noise)
    )
