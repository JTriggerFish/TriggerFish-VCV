"""FFT evaluation of the same truncated, reflected Gaussian power envelope."""

from functools import lru_cache
import numpy as np
from scipy.signal import fftconvolve


@lru_cache(maxsize=32)
def _kernel(sigma):
    radius = int(4 * sigma + 0.5)  # scipy.ndimage default truncate=4
    values = np.exp(-0.5 * (np.arange(-radius, radius + 1) / sigma) ** 2)
    values /= values.sum()
    return radius, values


def smoothed_power(samples, sigma_samples):
    """Equivalent to gaussian_filter1d(samples**2, sigma, mode='reflect')."""
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 1 or not samples.size or not np.isfinite(samples).all():
        raise ValueError("Expected finite nonempty mono audio")
    if not np.isfinite(sigma_samples) or sigma_samples <= 0:
        raise ValueError("Gaussian sigma must be positive")
    radius, kernel = _kernel(float(sigma_samples))
    padded = np.pad(samples**2, (radius, radius), mode="symmetric")
    return np.maximum(0, fftconvolve(padded, kernel, mode="valid"))
