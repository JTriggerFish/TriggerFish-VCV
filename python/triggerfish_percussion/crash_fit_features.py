"""Noise-aware perceptual features and weights for crash calibration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer
from .crash_fit_common import fixed_rate
from .descriptors import spectral_trajectories
from .erb import ErbFilterbank
from .transforms import StftConfig, stft


@dataclass(frozen=True)
class FeatureSet:
    values: np.ndarray
    names: tuple[str, ...]


REGIONS = (
    ("contact", 0.000, 0.015),
    ("bloom", 0.015, 0.250),
    ("early", 0.250, 1.500),
    ("tail-a", 1.500, 4.000),
    ("tail-b", 4.000, 8.000),
)

TRAJECTORY_TIMES_SECONDS = (
    0.025,
    0.05,
    0.1,
    0.2,
    0.4,
    0.8,
    1.5,
    2.5,
    4.0,
    6.0,
    8.0,
    10.0,
)


def perceptual_features(audio: AudioBuffer) -> FeatureSet:
    mono = audio.mono()
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    signal = mono.samples[onset:]
    if signal.size < 2048:
        signal = np.pad(signal, (0, 2048 - signal.size))
    transform = stft(signal, mono.sample_rate, StftConfig(2048, 512, 4096))
    bank = ErbFilterbank.create(transform.frequencies_hz, 80.0, 18000.0, 64)
    noise_power = _stationary_noise_power(transform.power, transform.times_seconds)
    smoothed_power = gaussian_filter1d(transform.power, 2.0, axis=1, mode="nearest")
    clean_power = np.maximum(
        smoothed_power - noise_power[:, np.newaxis],
        max(float(np.max(transform.power)) * 1.0e-14, np.finfo(float).tiny),
    )
    raw_band_power = bank.apply_power(transform.power)
    noise_band_power = bank.apply_power(noise_power[:, np.newaxis])[:, 0]
    band_floor = max(float(np.max(raw_band_power)) * 1.0e-12, np.finfo(float).tiny)
    band_power = np.maximum(
        gaussian_filter1d(raw_band_power, 2.0, axis=1, mode="nearest")
        - noise_band_power[:, np.newaxis],
        band_floor,
    )
    density_bins = (transform.frequencies_hz >= 300.0) & (
        transform.frequencies_hz <= 18000.0
    )
    trajectories = spectral_trajectories(
        transform.frequencies_hz[density_bins], clean_power[density_bins]
    )
    values: list[float] = []
    names: list[str] = []
    floor = max(float(np.max(band_power)) * 1.0e-10, np.finfo(float).tiny)
    for region, start, end in REGIONS:
        selected = (transform.times_seconds >= start) & (transform.times_seconds < end)
        if not np.any(selected):
            region_power = np.full(band_power.shape[0], floor)
        else:
            region_power = np.maximum(
                np.mean(raw_band_power[:, selected], axis=1) - noise_band_power,
                floor,
            )
        total = max(float(np.sum(region_power)), floor)
        values.append(10.0 * np.log10(total) / 6.0)
        names.append(f"level/{region}")
        shape_db = 10.0 * np.log10(np.maximum(region_power / total, 1.0e-10))
        values.extend((shape_db / 12.0).tolist())
        names.extend(f"erb/{region}/{frequency:.0f}" for frequency in bank.centers_hz)
        if np.any(selected):
            density_values = (
                np.log(max(float(np.median(trajectories.centroid_hz[selected])), 1.0))
                / 0.5,
                np.log(max(float(np.median(trajectories.flatness[selected])), 1.0e-8)),
                np.log(max(float(np.median(trajectories.crest[selected])), 1.0e-8)),
                float(np.median(trajectories.flux[selected])) / 0.05,
            )
        else:
            density_values = (0.0, 0.0, 0.0, 0.0)
        values.extend(density_values)
        names.extend(
            f"{descriptor}/{region}"
            for descriptor in ("centroid", "flatness", "crest", "flux")
        )
    trajectory_peak = max(
        float(np.max(np.sum(band_power, axis=0))), np.finfo(float).tiny
    )
    trajectory_floor = trajectory_peak * 1.0e-9
    for time_seconds in TRAJECTORY_TIMES_SECONDS:
        frame = int(np.argmin(np.abs(transform.times_seconds - time_seconds)))
        level_db = 10.0 * np.log10(
            np.maximum(band_power[:, frame], trajectory_floor) / trajectory_peak
        )
        values.extend((level_db / 12.0).tolist())
        names.extend(
            f"trajectory/{time_seconds:.3f}/{frequency:.0f}"
            for frequency in bank.centers_hz
        )
    persistent_bins = np.flatnonzero(
        (transform.frequencies_hz >= 100.0) & (transform.frequencies_hz <= 6000.0)
    )[::4]
    for region, start, end in (("early", 0.2, 1.5), ("tail", 1.5, 4.0)):
        selected = (transform.times_seconds >= start) & (transform.times_seconds < end)
        if np.any(selected):
            profile = gaussian_filter1d(
                np.mean(clean_power[:, selected], axis=1), 1.5, mode="nearest"
            )
            peak = max(float(np.max(profile[persistent_bins])), np.finfo(float).tiny)
            profile_db = 10.0 * np.log10(
                np.maximum(profile[persistent_bins], peak * 1.0e-8) / peak
            )
        else:
            profile_db = np.full(persistent_bins.size, -80.0)
        values.extend((profile_db / 12.0).tolist())
        names.extend(
            f"persistent/{region}/{transform.frequencies_hz[index]:.0f}"
            for index in persistent_bins
        )
    return FeatureSet(np.asarray(values), tuple(names))


def maximum_late_regrowth_db(audio: AudioBuffer) -> float:
    mono = fixed_rate(audio)
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    start = onset + round(2.0 * mono.sample_rate)
    window = round(0.25 * mono.sample_rate)
    levels = []
    for frame_start in range(start, mono.samples.size - window + 1, window):
        frame = mono.samples[frame_start : frame_start + window]
        levels.append(float(np.sqrt(np.mean(np.square(frame)))))
    if len(levels) < 2:
        return 0.0
    floor = max(max(levels) * 1.0e-9, np.finfo(float).tiny)
    levels_db = 20.0 * np.log10(np.maximum(levels, floor))
    prior_minimum = np.minimum.accumulate(levels_db)
    return float(max(0.0, np.max(levels_db[1:] - prior_minimum[:-1])))


def _stationary_noise_power(power: np.ndarray, times_seconds: np.ndarray) -> np.ndarray:
    if power.shape[1] == 0:
        return np.zeros(power.shape[0], dtype=np.float64)
    final_time = float(times_seconds[-1]) if times_seconds.size else 0.0
    selected = times_seconds >= max(0.0, final_time - 1.0)
    if np.count_nonzero(selected) < min(8, power.shape[1]):
        selected = np.zeros(power.shape[1], dtype=bool)
        selected[-min(8, power.shape[1]) :] = True
    return np.mean(power[:, selected], axis=1)
