"""Named contact and dense-response descriptors."""

from dataclasses import dataclass

import numpy as np
from scipy.signal import hilbert


@dataclass(frozen=True)
class ContactDescriptors:
    energy: float
    peak: float
    crest_factor: float
    spectral_centroid_hz: float
    spectral_bandwidth_hz: float


@dataclass(frozen=True)
class SpectralTrajectories:
    centroid_hz: np.ndarray
    bandwidth_hz: np.ndarray
    flatness: np.ndarray
    crest: np.ndarray
    flux: np.ndarray


def contact_descriptors(samples: np.ndarray, sample_rate: float) -> ContactDescriptors:
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or signal.size < 4 or not np.isfinite(signal).all():
        raise ValueError("contact descriptors require finite mono audio")
    envelope = np.abs(hilbert(signal))
    energy = float(np.sum(np.square(signal)))
    rms = float(np.sqrt(np.mean(np.square(signal))))
    peak = float(np.max(envelope))
    frequencies = np.fft.rfftfreq(signal.size, 1.0 / sample_rate)
    power = np.square(np.abs(np.fft.rfft(signal * np.hanning(signal.size))))
    centroid, bandwidth = _moments(frequencies, power[:, None])
    return ContactDescriptors(
        energy,
        peak,
        peak / max(rms, np.finfo(float).tiny),
        float(centroid[0]),
        float(bandwidth[0]),
    )


def spectral_trajectories(
    frequencies_hz: np.ndarray, power: np.ndarray
) -> SpectralTrajectories:
    frequencies = np.asarray(frequencies_hz, dtype=np.float64)
    values = np.asarray(power, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != frequencies.size or np.any(values < 0):
        raise ValueError("spectral trajectories require frequency-by-time power")
    centroid, bandwidth = _moments(frequencies, values)
    floor = np.maximum(values, np.finfo(float).tiny)
    flatness = np.exp(np.mean(np.log(floor), axis=0)) / np.mean(floor, axis=0)
    crest = np.max(values, axis=0) / np.maximum(
        np.mean(values, axis=0), np.finfo(float).tiny
    )
    normalized = values / np.maximum(
        np.sum(values, axis=0, keepdims=True), np.finfo(float).tiny
    )
    flux = np.sqrt(
        np.sum(
            np.square(np.diff(normalized, axis=1, prepend=normalized[:, :1])), axis=0
        )
    )
    return SpectralTrajectories(centroid, bandwidth, flatness, crest, flux)


def _moments(
    frequencies: np.ndarray, power: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    total = np.maximum(np.sum(power, axis=0), np.finfo(float).tiny)
    centroid = np.sum(frequencies[:, None] * power, axis=0) / total
    variance = (
        np.sum(np.square(frequencies[:, None] - centroid) * power, axis=0) / total
    )
    return centroid, np.sqrt(np.maximum(variance, 0.0))
