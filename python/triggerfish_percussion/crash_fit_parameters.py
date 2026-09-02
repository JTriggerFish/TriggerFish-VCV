"""Bounded scalar parameters and their causal crash fitting stages."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .crash_model import DENSE_GAIN_ENVELOPE_POINT_COUNT


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


BODY_DECAY_PARAMETERS = tuple(
    parameter(f"body_decay_seconds[{index}]", 0.02, 20.0) for index in range(5)
)
TURBULENCE_PARAMETERS = (
    *(parameter(f"turbulence_gain[{index}]", 0.0, 1.0) for index in range(3)),
    parameter("turbulence_persistence", 0.25, 4.0),
)


CAUSAL_STAGES = (
    CausalStage(
        "impact-contact",
        0.004,
        (
            parameter("contact_duration_scale", 0.35, 2.5),
            parameter("contact_pulse_gain", 0.0, 3.0),
            parameter("contact_chirp_gain", 0.0, 3.0),
            parameter("contact_chirp_frequency_scale", 0.5, 2.0),
            parameter("contact_noise_gain", 0.0, 4.0),
            parameter("contact_noise_duration_scale", 0.35, 3.0),
            parameter("contact_noise_tilt_db", -10.0, 10.0),
            parameter("contact_micro_gain", 0.0, 4.0),
            parameter("contact_micro_duration_scale", 0.35, 3.0),
            parameter("contact_micro_density_scale", 0.35, 2.5),
        ),
        requires_quality=False,
    ),
    CausalStage(
        "impact-balance",
        0.004,
        (
            parameter("direct_gain", 0.0, 2.0),
            parameter("sparse_gain", 0.0, 2.0),
            parameter("dense_gain", 0.0, 2.0),
            parameter("body_bypass_gain", 0.0, 0.5),
            parameter("strength_gamma", 0.35, 3.0),
            parameter("body_strength_gamma", 0.2, 2.5),
            parameter("dense_strength_gamma", 0.2, 2.5),
            parameter("dense_frequency_warp", 0.5, 2.0),
            parameter("dense_tilt_db_per_octave", -8.0, 5.0),
            parameter("high_cut_hz", 6000.0, 22000.0),
        ),
    ),
    CausalStage(
        "contact-tail",
        0.015,
        (
            parameter("contact_duration_scale", 0.35, 2.5),
            parameter("contact_pulse_gain", 0.0, 3.0),
            parameter("contact_chirp_gain", 0.0, 3.0),
            parameter("contact_chirp_frequency_scale", 0.5, 2.0),
            parameter("contact_noise_gain", 0.0, 4.0),
            parameter("contact_noise_duration_scale", 0.35, 3.0),
            parameter("contact_noise_tilt_db", -10.0, 10.0),
            parameter("contact_micro_gain", 0.0, 4.0),
            parameter("contact_micro_duration_scale", 0.35, 3.0),
            parameter("contact_micro_density_scale", 0.35, 2.5),
            parameter("dispersion_drive", 0.0, 8.0),
            parameter("dispersion_excursion_samples", 0.0, 16.0),
            parameter("dense_minimum_frequency_hz", 200.0, 1800.0),
            parameter("dense_maximum_frequency_hz", 12000.0, 22000.0),
            *(
                parameter(f"dense_gain_envelope_db[{index}]", -24.0, 24.0)
                for index in range(1, DENSE_GAIN_ENVELOPE_POINT_COUNT)
            ),
            parameter("dispersion_feedback", 0.7, 0.9995, 0.002),
            parameter("dispersion_low_decay_seconds", 0.1, 2.0),
            parameter("dispersion_middle_decay_seconds", 0.08, 1.5),
            parameter("dispersion_high_decay_seconds", 0.05, 1.0),
        ),
    ),
    CausalStage(
        "initial-decay",
        0.100,
        (
            parameter("sparse_gain", 0.0, 2.0),
            parameter("dense_gain", 0.0, 2.0),
            parameter("body_bypass_gain", 0.0, 0.5),
            parameter("dense_frequency_warp", 0.5, 2.0),
            parameter("dense_tilt_db_per_octave", -8.0, 5.0),
            *(
                parameter(f"dense_gain_envelope_db[{index}]", -24.0, 24.0)
                for index in range(1, DENSE_GAIN_ENVELOPE_POINT_COUNT)
            ),
            parameter("dense_minimum_frequency_hz", 200.0, 1800.0),
            parameter("dense_maximum_frequency_hz", 12000.0, 22000.0),
            parameter("dispersion_drive", 0.0, 8.0),
            parameter("dispersion_excursion_samples", 0.0, 16.0),
            parameter("dispersion_feedback", 0.7, 0.9995, 0.002),
            parameter("dispersion_low_decay_seconds", 0.1, 2.0),
            parameter("dispersion_middle_decay_seconds", 0.08, 1.5),
            parameter("dispersion_high_decay_seconds", 0.05, 1.0),
            *TURBULENCE_PARAMETERS,
            parameter("sparse_bloom_gain", 0.0, 1.0),
            *BODY_DECAY_PARAMETERS,
            parameter("colour_frequency_hz", 1000.0, 14000.0),
            parameter("colour_gain_db", -12.0, 12.0),
            parameter("high_cut_hz", 6000.0, 22000.0),
        ),
    ),
    CausalStage(
        "bloom",
        0.250,
        (
            parameter("dense_spacing_jitter", 0.0, 1.0),
            parameter("dense_decay_spread_octaves", 0.0, 1.25),
            parameter("dense_gain_spread_db", 0.0, 8.0),
        ),
    ),
    CausalStage(
        "early-body",
        1.500,
        (parameter("dense_velocity_loss_nepers_per_second", 0.0, 4.0, 0.1),),
    ),
    CausalStage(
        "tail-refinement",
        4.000,
        (*BODY_DECAY_PARAMETERS,),
    ),
)


ALL_CAUSAL_PARAMETERS = tuple(
    {item.name: item for stage in CAUSAL_STAGES for item in stage.parameters}.values()
)

ATTACK_PARAMETERS = tuple(
    {
        item.name: item
        for stage in CAUSAL_STAGES
        if stage.end_seconds <= 0.100
        for item in stage.parameters
    }.values()
)


def _selected_parameters(*names: str) -> tuple[FitParameter, ...]:
    available = {parameter.name: parameter for parameter in ALL_CAUSAL_PARAMETERS}
    return tuple(available[name] for name in names)


_SCREENED_CONTACT = _selected_parameters(
    "contact_duration_scale",
    "contact_noise_gain",
    "contact_noise_duration_scale",
    "contact_noise_tilt_db",
    "contact_pulse_gain",
    "contact_micro_gain",
)
_SCREENED_BALANCE = _selected_parameters(
    "direct_gain",
    "sparse_gain",
    "dense_gain",
    "strength_gamma",
    "body_strength_gamma",
    "dense_strength_gamma",
)
_SCREENED_COLOUR = _selected_parameters(
    "dense_tilt_db_per_octave",
    "dense_minimum_frequency_hz",
    "dense_maximum_frequency_hz",
    "sparse_bloom_gain",
    "colour_frequency_hz",
    "colour_gain_db",
    "high_cut_hz",
)
_SCREENED_RESIDUAL = _selected_parameters(
    "dispersion_drive",
    "dispersion_feedback",
    "dispersion_excursion_samples",
    "dense_frequency_warp",
    "contact_micro_density_scale",
    "turbulence_gain[0]",
    "turbulence_gain[1]",
    "turbulence_gain[2]",
    "turbulence_persistence",
)
_SCREENED_DECAY = _selected_parameters(
    *(f"body_decay_seconds[{index}]" for index in range(5)),
)
_SCREENED_JOINT = tuple(
    {
        parameter.name: parameter
        for group in (
            _SCREENED_CONTACT,
            _SCREENED_BALANCE,
            _SCREENED_COLOUR,
            _SCREENED_RESIDUAL,
            _SCREENED_DECAY,
        )
        for parameter in group
    }.values()
)

SCREENED_ATTACK_STAGES = (
    CausalStage(
        "screened-contact",
        0.100,
        _SCREENED_CONTACT,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "screened-balance",
        0.100,
        _SCREENED_BALANCE,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "screened-colour",
        0.100,
        _SCREENED_COLOUR,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "screened-residual",
        0.100,
        _SCREENED_RESIDUAL,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "screened-decay",
        0.100,
        _SCREENED_DECAY,
        requires_acceptance_gate=False,
    ),
    CausalStage("screened-joint-polish", 0.100, _SCREENED_JOINT),
)


_INITIAL_DECAY_BALANCE = _selected_parameters(
    "sparse_gain",
    "dense_gain",
    "body_bypass_gain",
    "sparse_bloom_gain",
)
_INITIAL_DECAY_COLOUR = _selected_parameters(
    "dense_frequency_warp",
    "dense_tilt_db_per_octave",
    *(
        f"dense_gain_envelope_db[{index}]"
        for index in range(1, DENSE_GAIN_ENVELOPE_POINT_COUNT)
    ),
    "dense_minimum_frequency_hz",
    "dense_maximum_frequency_hz",
    "colour_frequency_hz",
    "colour_gain_db",
    "high_cut_hz",
)
_INITIAL_DECAY_DISPERSION = _selected_parameters(
    "dispersion_drive",
    "dispersion_excursion_samples",
    "dispersion_feedback",
    "dispersion_low_decay_seconds",
    "dispersion_middle_decay_seconds",
    "dispersion_high_decay_seconds",
)
_INITIAL_DECAY_TURBULENCE = _selected_parameters(
    "turbulence_gain[0]",
    "turbulence_gain[1]",
    "turbulence_gain[2]",
    "turbulence_persistence",
)
_INITIAL_DECAY_LOSS = _selected_parameters(
    *(f"body_decay_seconds[{index}]" for index in range(5)),
)
_INITIAL_DECAY_JOINT = tuple(
    {
        item.name: item
        for group in (
            _INITIAL_DECAY_BALANCE,
            _INITIAL_DECAY_COLOUR,
            _INITIAL_DECAY_DISPERSION,
            _INITIAL_DECAY_TURBULENCE,
            _INITIAL_DECAY_LOSS,
        )
        for item in group
    }.values()
)

SCREENED_INITIAL_DECAY_STAGES = CAUSAL_STAGES[:3] + (
    CausalStage(
        "initial-decay-loss",
        0.100,
        _INITIAL_DECAY_LOSS,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "initial-decay-balance",
        0.100,
        _INITIAL_DECAY_BALANCE,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "initial-decay-colour",
        0.100,
        _INITIAL_DECAY_COLOUR,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "initial-decay-dispersion",
        0.100,
        _INITIAL_DECAY_DISPERSION,
        requires_acceptance_gate=False,
    ),
    CausalStage(
        "initial-decay-turbulence",
        0.100,
        _INITIAL_DECAY_TURBULENCE,
        requires_acceptance_gate=False,
    ),
    CausalStage("initial-decay-joint-polish", 0.100, _INITIAL_DECAY_JOINT),
)


SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS = frozenset(
    {
        "strength_gamma",
        "body_strength_gamma",
        "dense_strength_gamma",
        "dense_velocity_loss_nepers_per_second",
    }
)

SINGLE_HIT_ATTACK_PARAMETERS = tuple(
    parameter
    for parameter in ATTACK_PARAMETERS
    if parameter.name not in SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS
)


def single_hit_stages(
    stages: tuple[CausalStage, ...],
) -> tuple[CausalStage, ...]:
    """Remove response-curve controls that one reference hit cannot identify."""
    result = []
    for stage in stages:
        identifiable = tuple(
            parameter
            for parameter in stage.parameters
            if parameter.name not in SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS
        )
        if identifiable:
            result.append(replace(stage, parameters=identifiable))
    return tuple(result)
