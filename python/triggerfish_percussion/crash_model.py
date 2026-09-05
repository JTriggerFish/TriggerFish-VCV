"""Thin Python facade over the production C++ crash engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields

import numpy as np

BODY_DECAY_POINT_COUNT = 8
BODY_DECAY_INTERIOR_POINT_COUNT = BODY_DECAY_POINT_COUNT - 2


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
        15000.0,
        15000.0,
        15000.0,
        15000.0,
        15000.0,
        15000.0,
        15000.0,
        15000.0,
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
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
    field_turbulence_scale: tuple[float, ...] = (1.0,) * 32
    sparse_tune: float = 1.0
    body_decay_frequency_hz: tuple[float, ...] = (
        500.0,
        1500.0,
        5000.0,
        8000.0,
        12000.0,
        14000.0,
    )
    body_decay_seconds: tuple[float, ...] = (
        4.0,
        2.965782,
        2.372841,
        1.782984,
        1.585960,
        1.431738,
        1.330840,
        1.2,
    )
    body_decay_active: tuple[bool, ...] = (
        False,
        False,
        False,
        False,
        False,
        False,
    )
    body_tilt_db_per_octave: float = -1.0
    body_excitation_centre_hz: float = 1000.0
    bloom_rate_octaves_per_second: float = 2.0
    bloom_energy_acceleration: float = 0.7
    bloom_phase_diffusion: float = 0.7
    body_excitation_gain: float = 1.0
    field_gain: float = 1.0
    field_turbulence: float = 0.65
    field_turbulence_slope_per_octave: float = 0.0
    field_turbulence_centre_hz: float = 4000.0
    field_packet_spread_erb: float = 6.0
    field_satellite_density: float = 0.5
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
    output_gain: float = 1.0
    direct_radiation_enabled: bool = True
    direct_low_cut_hz: float = 40.0
    direct_colour_frequency_hz: float = 7200.0
    direct_colour_gain_db: float = 1.0
    direct_high_cut_hz: float = 20000.0
    body_radiation_enabled: bool = True
    body_low_cut_hz: float = 40.0
    body_colour_frequency_hz: float = 7200.0
    body_colour_gain_db: float = 0.5
    body_high_cut_hz: float = 19000.0
    velocity_brightness_db_per_octave: float = 4.0

    def __post_init__(self):
        frequencies = tuple(float(value) for value in self.body_decay_frequency_hz)
        seconds = tuple(float(value) for value in self.body_decay_seconds)
        active = tuple(bool(value) for value in self.body_decay_active)
        if (
            len(frequencies) != BODY_DECAY_INTERIOR_POINT_COUNT
            or len(active) != BODY_DECAY_INTERIOR_POINT_COUNT
            or len(seconds) != BODY_DECAY_POINT_COUNT
        ):
            raise ValueError("body decay curve needs six interior and eight T60 points")
        object.__setattr__(self, "body_decay_frequency_hz", frequencies)
        object.__setattr__(self, "body_decay_seconds", seconds)
        object.__setattr__(self, "body_decay_active", active)

    def native(self):
        import _triggerfish_dsp as native

        parameters = native.CrashCymbalFitParameters()
        for name, value in asdict(self).items():
            setattr(parameters, name, value)
        return parameters

    @classmethod
    def from_mapping(cls, values: dict[str, object]) -> "CrashFit":
        """Load the current fit schema without implicit legacy translation."""
        current = dict(values)
        known = {item.name for item in fields(cls)}
        unknown = sorted(set(current) - known)
        if unknown:
            raise ValueError(f"unknown crash fit fields: {', '.join(unknown)}")
        return cls(**current)


@dataclass(frozen=True)
class CrashEvent:
    onset_seconds: float
    strength: float
    location: float = 1.0
    hardness: float = 0.65
    seed: int = 1
    implement: float = 0.75
    contact_spread: float = 0.2


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
    """Return direct contact, cascade transfer, modal body, and output taps."""
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
