"""Bounded parameters for the active unified metallic-body fitter."""

from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class FitParameter:
    name: str
    lower: float
    upper: float
    relative_step: float = 0.01


@dataclass(frozen=True)
class CausalStage:
    name: str
    end_seconds: float
    parameters: tuple[FitParameter, ...]
    requires_quality: bool = True
    requires_acceptance_gate: bool = True
    maximum_envelope_rmse_db: float = 6.0
    maximum_envelope_absolute_db: float = 6.0
    maximum_spectral_rmse_db: float = 10.0
    maximum_spectral_p95_db: float = 18.0
    maximum_fine_spectrum_rmse_db: float = 4.0
    maximum_centroid_rmse_octaves: float = 0.35
    maximum_rolloff_rmse_octaves: float = 0.35
    maximum_flatness_rmse_db: float = 2.5
    maximum_crest_rmse_db: float = 2.5
    maximum_ridge_ratio_error: float = 0.08


def parameter(
    name: str, lower: float, upper: float, step: float = 0.01
) -> FitParameter:
    return FitParameter(name, lower, upper, step)


def fit_parameter_value(fit, name: str) -> float:
    """Read a scalar or one indexed tuple field from a frozen fit object."""
    field, index = _parameter_address(name)
    value = getattr(fit, field)
    return float(value if index is None else value[index])


def replace_fit_parameters(fit, values: dict[str, float]):
    """Replace scalar and indexed tuple fields in one dataclass operation."""
    changes = {}
    for name, value in values.items():
        field, index = _parameter_address(name)
        if index is None:
            changes[field] = float(value)
            continue
        sequence = list(changes.get(field, getattr(fit, field)))
        sequence[index] = float(value)
        changes[field] = tuple(sequence)
    return replace(fit, **changes)


def _parameter_address(name: str) -> tuple[str, int | None]:
    if "[" not in name:
        return name, None
    field, suffix = name.split("[", 1)
    if not suffix.endswith("]"):
        raise ValueError(f"invalid indexed fit parameter: {name}")
    return field, int(suffix[:-1])


# Ordinary fitting owns only the fixed DC and Nyquist values. Interior slots
# remain inactive until an explicit last-step refinement is justified.
BODY_DECAY_SLOTS = (0, 7)
BODY_DECAY_PARAMETERS = tuple(
    parameter(f"body_decay_seconds[{index}]", 0.02, 20.0) for index in BODY_DECAY_SLOTS
)


UNIFIED_CAUSAL_STAGES = (
    CausalStage(
        "impact-contact",
        0.004,
        (
            parameter("contact_duration_scale", 0.25, 4.0),
            parameter("contact_pulse_gain", 0.0, 4.0),
            parameter("contact_chirp_gain", 0.0, 4.0),
            parameter("contact_chirp_frequency_scale", 0.05, 4.0),
            parameter("contact_noise_gain", 0.0, 4.0),
            parameter("contact_noise_duration_scale", 0.25, 4.0),
            parameter("contact_noise_tilt_db", -18.0, 18.0),
            parameter("contact_micro_gain", 0.0, 4.0),
            parameter("contact_micro_duration_scale", 0.25, 4.0),
            parameter("contact_micro_density_scale", 0.25, 4.0),
        ),
        requires_quality=False,
    ),
    CausalStage(
        "unified-impact-balance",
        0.015,
        (
            parameter("direct_gain", 0.0, 2.0),
            parameter("body_excitation_gain", 0.001, 4.0),
            parameter("field_gain", 0.0, 3.0),
            parameter("direct_colour_frequency_hz", 100.0, 18_000.0),
            parameter("direct_colour_gain_db", -18.0, 18.0),
            parameter("direct_high_cut_hz", 1_000.0, 22_000.0),
        ),
    ),
    CausalStage(
        "unified-initial-body",
        0.100,
        (
            parameter("field_gain", 0.0, 3.0),
            parameter("body_excitation_gain", 0.001, 4.0),
            parameter("body_tilt_db_per_octave", -24.0, 24.0),
            parameter("body_tilt_centre_hz", 40.0, 15_000.0),
            parameter("field_turbulence", 0.0, 1.0),
            parameter("field_turbulence_slope_per_octave", -1.0, 1.0),
            parameter("field_turbulence_centre_hz", 40.0, 15_000.0),
            parameter("field_packet_spread_erb", 0.0, 12.0),
            parameter("field_phase_bandwidth_erb", 0.0, 4.0),
            parameter("field_exchange", 0.0, 1.0),
            parameter("bloom_rate_octaves_per_second", 0.0, 16.0),
            parameter("bloom_energy_dependence", 0.0, 1.0),
            parameter("bloom_phase_diffusion", 0.0, 1.0),
            parameter("body_colour_frequency_hz", 100.0, 18_000.0),
            parameter("body_colour_gain_db", -18.0, 18.0),
            parameter("body_high_cut_hz", 1_000.0, 22_000.0),
        ),
    ),
    CausalStage(
        "unified-bloom",
        0.600,
        (
            parameter("bloom_rate_octaves_per_second", 0.0, 16.0),
            parameter("bloom_energy_dependence", 0.0, 1.0),
            parameter("bloom_phase_diffusion", 0.0, 1.0),
            parameter("field_turbulence", 0.0, 1.0),
            parameter("field_turbulence_slope_per_octave", -1.0, 1.0),
            parameter("field_turbulence_centre_hz", 40.0, 15_000.0),
            parameter("field_packet_spread_erb", 0.0, 12.0),
            parameter("field_phase_bandwidth_erb", 0.0, 4.0),
            parameter("field_exchange", 0.0, 1.0),
        ),
    ),
    CausalStage("unified-tail", 4.0, BODY_DECAY_PARAMETERS),
)

# Historical exported names resolve to the one renderer topology. There is no
# second parameter universe containing disconnected sparse/dense branches.
CAUSAL_STAGES = UNIFIED_CAUSAL_STAGES
ALL_CAUSAL_PARAMETERS = tuple(
    {
        item.name: item for stage in UNIFIED_CAUSAL_STAGES for item in stage.parameters
    }.values()
)
ATTACK_PARAMETERS = tuple(
    {
        item.name: item
        for stage in UNIFIED_CAUSAL_STAGES
        if stage.end_seconds <= 0.100
        for item in stage.parameters
    }.values()
)
SCREENED_ATTACK_STAGES = UNIFIED_CAUSAL_STAGES[:3]
SCREENED_INITIAL_DECAY_STAGES = UNIFIED_CAUSAL_STAGES[:4]
SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS = frozenset()
SINGLE_HIT_ATTACK_PARAMETERS = ATTACK_PARAMETERS


def single_hit_stages(
    stages: tuple[CausalStage, ...],
) -> tuple[CausalStage, ...]:
    """Return the active stages; no hidden response-curve controls remain."""
    return tuple(replace(stage, parameters=tuple(stage.parameters)) for stage in stages)
