"""Thin Python facade over the production C++ crash engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

DENSE_GAIN_ENVELOPE_POINT_COUNT = 33
BODY_DECAY_POINT_COUNT = 8


@dataclass(frozen=True)
class CrashFit:
    sparse_frequency_hz: tuple[float, ...] = (
        421.0,
        522.0,
        689.0,
        1094.0,
        1475.0,
        2009.0,
        2138.0,
        2573.0,
        2753.0,
        3589.0,
        4428.0,
        5707.0,
        6500.0,
        7350.0,
        8250.0,
        9200.0,
        10250.0,
        11350.0,
        12000.0,
        12750.0,
        13500.0,
        14100.0,
        14600.0,
        15000.0,
    )
    sparse_decay_ratio: tuple[float, ...] = (
        0.7,
        0.7,
        0.7,
        1.25,
        1.0,
        1.0,
        0.95,
        0.85,
        0.8,
        0.7,
        0.7,
        0.7,
        0.68,
        0.65,
        0.62,
        0.6,
        0.58,
        0.56,
        0.55,
        0.54,
        0.53,
        0.52,
        0.51,
        0.5,
    )
    sparse_amplitude: tuple[float, ...] = (
        0.35,
        0.15,
        0.15,
        0.55,
        0.25,
        0.7,
        0.65,
        0.7,
        0.55,
        0.5,
        0.4,
        0.25,
        0.32,
        0.38,
        0.44,
        0.5,
        0.56,
        0.6,
        0.62,
        0.6,
        0.54,
        0.46,
        0.35,
        0.22,
    )
    sparse_phase_radians: tuple[float, ...] = (0.0,) * 24
    field_turbulence_scale: tuple[float, ...] = (1.0,) * 24
    sparse_tune: float = 1.0
    body_decay_frequency_hz: tuple[float, ...] = (
        0.0,
        500.0,
        1500.0,
        5000.0,
        8000.0,
        12000.0,
        16000.0,
        24000.0,
    )
    body_decay_seconds: tuple[float, ...] = (4.0, 4.0, 3.8, 2.3, 1.8, 1.5, 1.3, 1.2)
    body_decay_active: tuple[bool, ...] = (
        True,
        True,
        True,
        True,
        False,
        False,
        False,
        True,
    )
    dense_minimum_frequency_hz: float = 180.0
    dense_maximum_frequency_hz: float = 18000.0
    dense_frequency_warp: float = 1.0
    dense_spacing_jitter: float = 0.82
    dense_mode_density: float = 1.0
    dense_decay_spread_octaves: float = 0.15
    dense_tilt_db_per_octave: float = -1.0
    dense_gain_envelope_db: tuple[float, ...] = (
        4.0,
        1.375,
        -1.25,
        -3.875,
        -6.5,
        -7.25,
        -5.5,
        -3.75,
        -2.0,
        -0.25,
        1.125,
        2.4375,
        3.75,
        5.0625,
        6.125,
        6.5625,
        7.0,
        7.4375,
        7.875,
        7.6875,
        7.25,
        6.8125,
        6.375,
        5.9375,
        5.5,
        5.0625,
        4.625,
        4.1875,
        3.75,
        3.3125,
        2.875,
        2.4375,
        2.0,
    )
    dense_gain_spread_db: float = 2.0
    dense_mode_seed: int = 0x43524153
    turbulence_frequency_hz: tuple[float, ...] = (350.0, 2200.0, 10000.0)
    turbulence_gain: tuple[float, ...] = (0.14, 0.1, 0.06)
    turbulence_persistence: float = 1.0
    dispersion_feedback: float = 0.9965
    dispersion_drive: float = 2.8
    dispersion_excursion_samples: float = 2.4
    dispersion_low_decay_seconds: float = 0.9
    dispersion_middle_decay_seconds: float = 0.65
    dispersion_high_decay_seconds: float = 0.42
    dispersion_diffusion: float = 1.0
    bloom_body_gain: float = 1.0
    field_gain: float = 0.73824115
    field_turbulence: float = 0.65
    field_packet_spread_erb: float = 6.0
    field_phase_bandwidth_erb: float = 1.0
    field_exchange: float = 0.35
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
    direct_radiation_enabled: bool = True
    direct_low_cut_hz: float = 40.0
    direct_low_cut_q: float = 0.707
    direct_colour_frequency_hz: float = 7200.0
    direct_colour_gain_db: float = 1.0
    direct_colour_q: float = 0.8
    direct_high_cut_hz: float = 20000.0
    direct_high_cut_q: float = 0.707
    sparse_radiation_enabled: bool = True
    sparse_low_cut_hz: float = 40.0
    sparse_low_cut_q: float = 0.707
    colour_frequency_hz: float = 5200.0
    colour_gain_db: float = 1.5
    sparse_colour_q: float = 0.8
    high_cut_hz: float = 19000.0
    sparse_high_cut_q: float = 0.707
    dense_radiation_enabled: bool = True
    dense_low_cut_hz: float = 40.0
    dense_low_cut_q: float = 0.707
    dense_colour_frequency_hz: float = 7200.0
    dense_colour_gain_db: float = 0.5
    dense_colour_q: float = 0.8
    dense_high_cut_hz: float = 19000.0
    dense_high_cut_q: float = 0.707
    strength_gamma: float = 1.0
    body_strength_gamma: float = 1.0

    def __post_init__(self):
        object.__setattr__(
            self,
            "dense_gain_envelope_db",
            _resample_envelope(
                self.dense_gain_envelope_db, DENSE_GAIN_ENVELOPE_POINT_COUNT
            ),
        )
        frequencies, seconds, active = _upgrade_decay_curve(
            self.body_decay_frequency_hz,
            self.body_decay_seconds,
            self.body_decay_active,
        )
        object.__setattr__(self, "body_decay_frequency_hz", frequencies)
        object.__setattr__(self, "body_decay_seconds", seconds)
        object.__setattr__(self, "body_decay_active", active)

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
    implement: float = 0.75
    contact_spread: float = 0.2


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


def _upgrade_decay_curve(frequencies, seconds, active):
    """Upgrade the former five-knot curve to eight fixed-capacity slots."""
    frequencies = tuple(float(value) for value in frequencies)
    seconds = tuple(float(value) for value in seconds)
    active = tuple(bool(value) for value in active)
    if len(frequencies) == len(seconds) == len(active) == BODY_DECAY_POINT_COUNT:
        return frequencies, seconds, active
    if len(frequencies) == len(seconds) == 5:
        return (
            (0.0, *frequencies, 20000.0, 24000.0),
            (seconds[0], *seconds, seconds[-1], seconds[-1]),
            (True, True, True, True, True, True, False, True),
        )
    if len(frequencies) == BODY_DECAY_POINT_COUNT and len(seconds) == 5:
        seconds = (
            seconds[0],
            seconds[1],
            seconds[2],
            seconds[3],
            1.0,
            1.0,
            1.0,
            seconds[4],
        )
        return frequencies, seconds, (True, True, True, True, False, False, False, True)
    raise ValueError("body decay curve needs five legacy or eight current points")


def render_crash(
    fit: CrashFit,
    seconds: float,
    sample_rate: int = 48000,
    strength: float = 0.8,
    location: float = 1.0,
    hardness: float = 0.65,
    seed: int = 1,
    implement: float = 0.75,
    contact_spread: float = 0.2,
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
            implement,
            contact_spread,
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
    implement: float = 0.75,
    contact_spread: float = 0.2,
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
            implement,
            contact_spread,
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
            np.asarray([event.implement for event in ordered], dtype=np.float32),
            np.asarray([event.contact_spread for event in ordered], dtype=np.float32),
            onsets,
            np.asarray([event.seed for event in ordered], dtype=np.uint32),
            fit.native(),
        ),
        dtype=np.float64,
    )
