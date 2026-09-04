"""Shared crash-fitting cells, rendering helpers, and analysis windows."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer, resample_audio
from .crash_model import CrashFit, render_crash
from .erb import ErbFilterbank
from .transforms import StftConfig, stft


@dataclass(frozen=True)
class CrashFitCell:
    label: str
    reference: AudioBuffer
    strength: float
    location: float = 1.0
    hardness: float = 0.65
    seed: int = 1


def fixed_rate(audio: AudioBuffer) -> AudioBuffer:
    return resample_audio(audio.mono(), 48000)


def render_cell(cell: CrashFitCell, fit: CrashFit, seconds: float) -> AudioBuffer:
    return AudioBuffer(
        render_crash(
            fit,
            seconds,
            48000,
            cell.strength,
            cell.location,
            cell.hardness,
            cell.seed,
        ),
        48000,
    )


def render_sparse(cell: CrashFitCell, fit: CrashFit, seconds: float) -> AudioBuffer:
    return render_cell(cell, replace(fit, direct_gain=0.0), seconds)


def body_rms(audio: AudioBuffer) -> float:
    mono = fixed_rate(audio)
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    return window_rms(mono, onset, 0.25, 1.5)


def window_rms(
    audio: AudioBuffer, onset: int, start_seconds: float, end_seconds: float
) -> float:
    start = min(audio.samples.size, onset + round(start_seconds * audio.sample_rate))
    end = min(audio.samples.size, onset + round(end_seconds * audio.sample_rate))
    if end <= start:
        return 0.0
    return float(np.sqrt(np.mean(np.square(audio.samples[start:end]))))


def modal_power_features(audio: AudioBuffer, onset: int) -> np.ndarray:
    """Return phase-insensitive, time-resolved modal-band energies."""
    start = onset + round(0.22 * audio.sample_rate)
    sample_count = round((2.2 - 0.22) * audio.sample_rate)
    window = np.zeros(sample_count, dtype=np.float64)
    available = max(0, min(sample_count, audio.samples.size - start))
    if available:
        window[:available] = audio.samples[start : start + available]
    transform = stft(window, audio.sample_rate, StftConfig(4096, 2048, 8192))
    bank = ErbFilterbank.create(
        transform.frequencies_hz,
        40.0,
        min(15_000.0, 0.48 * audio.sample_rate),
        96,
    )
    band_power = bank.apply_power(transform.power)
    # Averaging adjacent frames suppresses finite-window cross terms without
    # discarding the decay trajectory that identifies modal energy.
    frame_group = 4
    groups = []
    for first in range(0, band_power.shape[1], frame_group):
        groups.append(np.mean(band_power[:, first : first + frame_group], axis=1))
    return np.concatenate(groups)
