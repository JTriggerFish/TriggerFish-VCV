"""Thin Python facade over the production C++ crash engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class CrashFit:
    resonance_tune: float = 1.0
    low_decay_scale: float = 1.0
    middle_decay_scale: float = 1.0
    high_decay_scale: float = 1.0
    resonator_coupling: float = 0.38
    resonator_shift_scale: float = 1.0
    dispersion_feedback: float = 0.93
    dispersion_drive: float = 2.8
    dispersion_excursion_samples: float = 2.4
    dispersion_low_decay_seconds: float = 0.9
    dispersion_middle_decay_seconds: float = 0.65
    dispersion_high_decay_seconds: float = 0.42
    direct_gain: float = 0.18
    body_gain: float = 0.12
    body_bypass_gain: float = 0.06
    output_gain: float = 1.0
    colour_frequency_hz: float = 5200.0
    colour_gain_db: float = 1.5
    high_cut_hz: float = 19000.0
    strength_gamma: float = 1.15

    def native(self):
        import _triggerfish_dsp as native

        parameters = native.CrashCymbalFitParameters()
        for name, value in asdict(self).items():
            setattr(parameters, name, value)
        return parameters


@dataclass(frozen=True)
class CrashEvent:
    onset_seconds: float
    strength: float
    location: float = 1.0
    hardness: float = 0.65
    seed: int = 1


def render_crash(
    fit: CrashFit,
    seconds: float,
    sample_rate: int = 48000,
    strength: float = 0.8,
    location: float = 1.0,
    hardness: float = 0.65,
    seed: int = 1,
) -> np.ndarray:
    if seconds <= 0 or sample_rate < 1:
        raise ValueError("crash render duration and sample rate must be positive")
    import _triggerfish_dsp as native

    return np.asarray(
        native.render_crash(
            round(seconds * sample_rate),
            sample_rate,
            strength,
            location,
            hardness,
            seed,
            fit.native(),
        ),
        dtype=np.float64,
    )


def render_crash_sequence(
    fit: CrashFit,
    seconds: float,
    events: tuple[CrashEvent, ...],
    sample_rate: int = 48000,
) -> np.ndarray:
    if seconds <= 0 or sample_rate < 1:
        raise ValueError("crash render duration and sample rate must be positive")
    ordered = tuple(sorted(events, key=lambda event: event.onset_seconds))
    onsets = np.asarray(
        [round(event.onset_seconds * sample_rate) for event in ordered],
        dtype=np.intp,
    )
    if np.any(onsets < 0) or np.any(onsets >= round(seconds * sample_rate)):
        raise ValueError("crash event onsets must lie inside the render")
    import _triggerfish_dsp as native

    return np.asarray(
        native.render_crash_sequence(
            round(seconds * sample_rate),
            sample_rate,
            np.asarray([event.strength for event in ordered], dtype=np.float32),
            np.asarray([event.location for event in ordered], dtype=np.float32),
            np.asarray([event.hardness for event in ordered], dtype=np.float32),
            onsets,
            np.asarray([event.seed for event in ordered], dtype=np.uint32),
            fit.native(),
        ),
        dtype=np.float64,
    )
