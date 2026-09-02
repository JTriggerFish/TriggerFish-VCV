"""Thin Python facade over the production C++ crash engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

DENSE_GAIN_ENVELOPE_POINT_COUNT = 33


@dataclass(frozen=True)
class CrashFit:
    sparse_frequency_hz: tuple[float, ...] = (
        522.0,
        689.0,
        1094.0,
        1475.0,
        2009.0,
        2138.0,
        2573.0,
        2753.0,
        3589.0,
        3923.0,
        4428.0,
        5707.0,
    )
    sparse_decay_seconds: tuple[float, ...] = (
        5.5,
        6.5,
        4.8,
        3.8,
        3.2,
        3.0,
        2.5,
        2.2,
        1.6,
        1.35,
        1.05,
        0.75,
    )
    sparse_amplitude: tuple[float, ...] = (
        0.25,
        0.32,
        0.55,
        0.25,
        0.8,
        0.65,
        0.75,
        0.62,
        0.55,
        0.48,
        0.4,
        0.25,
    )
    sparse_phase_radians: tuple[float, ...] = (0.0,) * 12
    sparse_tune: float = 1.0
    sparse_decay_scale: float = 1.0
    dense_minimum_frequency_hz: float = 700.0
    dense_maximum_frequency_hz: float = 18000.0
    dense_frequency_warp: float = 1.0
    dense_spacing_jitter: float = 0.82
    dense_low_decay_seconds: float = 3.2
    dense_high_decay_seconds: float = 0.22
    dense_decay_curve: float = 0.75
    dense_decay_envelope_octaves: tuple[float, ...] = (0.0,) * 6
    dense_decay_spread_octaves: float = 0.4
    dense_tilt_db_per_octave: float = -1.0
    dense_gain_envelope_db: tuple[float, ...] = (0.0,) * DENSE_GAIN_ENVELOPE_POINT_COUNT
    dense_gain_spread_db: float = 4.5
    dense_mode_seed: int = 0x43524153
    turbulence_low_gain: float = 0.0
    turbulence_middle_gain: float = 0.0
    turbulence_high_gain: float = 0.0
    turbulence_low_decay_seconds: float = 0.8
    turbulence_middle_decay_seconds: float = 0.55
    turbulence_high_decay_seconds: float = 0.3
    dispersion_feedback: float = 0.9965
    dispersion_drive: float = 2.8
    dispersion_excursion_samples: float = 2.4
    dispersion_low_decay_seconds: float = 0.9
    dispersion_middle_decay_seconds: float = 0.65
    dispersion_high_decay_seconds: float = 0.42
    contact_duration_scale: float = 1.0
    contact_pulse_gain: float = 1.0
    contact_chirp_gain: float = 1.0
    contact_chirp_frequency_scale: float = 1.0
    contact_noise_gain: float = 1.0
    contact_noise_duration_scale: float = 1.0
    contact_noise_tilt_db: float = 0.0
    contact_micro_gain: float = 1.0
    contact_micro_duration_scale: float = 1.0
    contact_micro_density_scale: float = 1.0
    direct_gain: float = 0.18
    sparse_gain: float = 0.35
    dense_gain: float = 0.65
    sparse_bloom_gain: float = 0.0
    body_bypass_gain: float = 0.06
    output_gain: float = 1.0
    colour_frequency_hz: float = 5200.0
    colour_gain_db: float = 1.5
    high_cut_hz: float = 19000.0
    strength_gamma: float = 1.15
    body_strength_gamma: float = 0.8
    dense_strength_gamma: float = 0.8
    dense_velocity_loss_nepers_per_second: float = 0.0

    def __post_init__(self):
        object.__setattr__(
            self,
            "dense_gain_envelope_db",
            _resample_envelope(
                self.dense_gain_envelope_db, DENSE_GAIN_ENVELOPE_POINT_COUNT
            ),
        )
        object.__setattr__(
            self,
            "dense_decay_envelope_octaves",
            _resample_envelope(self.dense_decay_envelope_octaves, 6),
        )

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


def _resample_envelope(values, point_count: int) -> tuple[float, ...]:
    """Migrate older smooth-envelope grids without changing their endpoints."""
    source = np.asarray(values, dtype=np.float64)
    if source.ndim != 1 or source.size < 2 or not np.isfinite(source).all():
        raise ValueError("crash envelope needs at least two finite points")
    if source.size == point_count:
        return tuple(float(value) for value in source)
    positions = np.linspace(0.0, 1.0, source.size)
    target = np.linspace(0.0, 1.0, point_count)
    return tuple(float(value) for value in np.interp(target, positions, source))


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


def render_crash_components(
    fit: CrashFit,
    seconds: float,
    sample_rate: int = 48000,
    strength: float = 0.8,
    location: float = 1.0,
    hardness: float = 0.65,
    seed: int = 1,
) -> np.ndarray:
    """Return direct, dispersion, sparse, dense-residual, and output taps."""
    if seconds <= 0 or sample_rate < 1:
        raise ValueError("crash render duration and sample rate must be positive")
    import _triggerfish_dsp as native

    return np.asarray(
        native.render_crash_components(
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
