"""High-resolution damped-mode estimation and one-to-one matching."""

from dataclasses import dataclass

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.linalg import svd
from scipy.optimize import linear_sum_assignment


@dataclass(frozen=True)
class DampedMode:
    frequency_hz: float
    decay_seconds: float
    amplitude: float
    phase_radians: float


@dataclass(frozen=True)
class ModeMatch:
    reference_index: int
    candidate_index: int
    frequency_error_cents: float
    decay_error_log: float
    amplitude_error_db: float


@dataclass(frozen=True)
class ModeMatching:
    matches: tuple[ModeMatch, ...]
    missing_reference: tuple[int, ...]
    extra_candidate: tuple[int, ...]


@dataclass(frozen=True)
class ModelOrderEvidence:
    selected_mode_count: int
    candidate_mode_counts: np.ndarray
    mdl_values: np.ndarray
    singular_values: np.ndarray


def estimate_damped_modes(
    samples: np.ndarray,
    sample_rate: float,
    mode_count: int,
    pencil_samples: int = 512,
    maximum_analysis_samples: int = 8192,
) -> tuple[DampedMode, ...]:
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or signal.size < 16 or not np.isfinite(signal).all():
        raise ValueError("ESPRIT requires finite non-empty mono audio")
    if mode_count < 1 or sample_rate <= 0:
        raise ValueError("ESPRIT mode count and sample rate must be positive")
    signal = signal[:maximum_analysis_samples]
    pencil = min(max(2 * mode_count + 2, pencil_samples), signal.size // 2)
    hankel = sliding_window_view(signal, pencil).T
    left_vectors, _, _ = svd(hankel, full_matrices=False, check_finite=False)
    order = min(2 * mode_count, left_vectors.shape[1] - 1)
    subspace = left_vectors[:, :order]
    transition = np.linalg.lstsq(subspace[:-1], subspace[1:], rcond=None)[0]
    poles = np.linalg.eigvals(transition)
    coefficients = _fit_coefficients(signal, poles)
    modes = _positive_modes(poles, coefficients, sample_rate)
    return tuple(sorted(modes, key=lambda mode: mode.frequency_hz)[:mode_count])


def estimate_mode_count(
    samples: np.ndarray,
    maximum_mode_count: int,
    pencil_samples: int = 512,
    maximum_analysis_samples: int = 8192,
) -> ModelOrderEvidence:
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or signal.size < 16 or not np.isfinite(signal).all():
        raise ValueError("model-order evidence requires finite mono audio")
    if maximum_mode_count < 1:
        raise ValueError("maximum mode count must be positive")
    signal = signal[:maximum_analysis_samples]
    pencil = min(max(2 * maximum_mode_count + 2, pencil_samples), signal.size // 2)
    hankel = sliding_window_view(signal, pencil).T
    singular_values = svd(
        hankel, full_matrices=False, compute_uv=False, check_finite=False
    )
    component_orders = 2 * np.arange(maximum_mode_count + 1)
    component_orders = component_orders[component_orders < singular_values.size]
    mdl = np.array(
        [
            _mdl_value(singular_values, int(order), hankel.shape[1])
            for order in component_orders
        ]
    )
    selected = int(component_orders[int(np.argmin(mdl))] // 2)
    return ModelOrderEvidence(selected, component_orders // 2, mdl, singular_values)


def match_modes(
    reference: tuple[DampedMode, ...],
    candidate: tuple[DampedMode, ...],
    frequency_scale_cents: float = 50.0,
    decay_scale_log: float = 0.35,
    amplitude_scale_db: float = 6.0,
    maximum_frequency_error_cents: float = 300.0,
    maximum_cost: float = 8.0,
) -> ModeMatching:
    if (
        frequency_scale_cents <= 0
        or decay_scale_log <= 0
        or amplitude_scale_db <= 0
        or maximum_frequency_error_cents <= 0
        or maximum_cost <= 0
    ):
        raise ValueError("mode-matching scales and gates must be positive")
    if not reference or not candidate:
        return ModeMatching(
            (), tuple(range(len(reference))), tuple(range(len(candidate)))
        )
    costs = np.empty((len(reference), len(candidate)), dtype=np.float64)
    details: dict[tuple[int, int], tuple[float, float, float]] = {}
    for row, ref in enumerate(reference):
        for column, current in enumerate(candidate):
            cents = 1200.0 * np.log2(current.frequency_hz / ref.frequency_hz)
            decay = np.log(current.decay_seconds / ref.decay_seconds)
            amplitude = 20.0 * np.log10(
                max(current.amplitude, 1.0e-12) / max(ref.amplitude, 1.0e-12)
            )
            details[row, column] = (cents, decay, amplitude)
            costs[row, column] = np.sqrt(
                (cents / frequency_scale_cents) ** 2
                + (decay / decay_scale_log) ** 2
                + (amplitude / amplitude_scale_db) ** 2
            )
    rows, columns = linear_sum_assignment(costs)
    accepted = [
        (row, column)
        for row, column in zip(rows.tolist(), columns.tolist())
        if abs(details[row, column][0]) <= maximum_frequency_error_cents
        and costs[row, column] <= maximum_cost
    ]
    accepted_rows = {row for row, _ in accepted}
    accepted_columns = {column for _, column in accepted}
    matches = tuple(
        ModeMatch(row, column, *details[row, column]) for row, column in accepted
    )
    return ModeMatching(
        matches,
        tuple(sorted(set(range(len(reference))) - accepted_rows)),
        tuple(sorted(set(range(len(candidate))) - accepted_columns)),
    )


def resynthesize_modes(
    modes: tuple[DampedMode, ...], sample_rate: float, sample_count: int
) -> np.ndarray:
    if sample_rate <= 0 or sample_count < 0:
        raise ValueError("modal resynthesis dimensions must be non-negative")
    time = np.arange(sample_count, dtype=np.float64) / sample_rate
    output = np.zeros(sample_count, dtype=np.float64)
    for mode in modes:
        if mode.frequency_hz <= 0 or mode.decay_seconds <= 0 or mode.amplitude < 0:
            raise ValueError("modal parameters must describe a decaying component")
        output += (
            mode.amplitude
            * np.exp(-time / mode.decay_seconds)
            * np.cos(2.0 * np.pi * mode.frequency_hz * time + mode.phase_radians)
        )
    return output


def refit_mode_amplitudes(
    samples: np.ndarray,
    sample_rate: float,
    modes: tuple[DampedMode, ...],
) -> tuple[DampedMode, ...]:
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or not np.isfinite(signal).all() or sample_rate <= 0:
        raise ValueError("modal refit requires finite mono audio")
    if not modes:
        return ()
    time = np.arange(signal.size, dtype=np.float64) / sample_rate
    columns: list[np.ndarray] = []
    for mode in modes:
        envelope = np.exp(-time / mode.decay_seconds)
        phase = 2.0 * np.pi * mode.frequency_hz * time
        columns.extend((envelope * np.cos(phase), envelope * np.sin(phase)))
    coefficients = np.linalg.lstsq(np.column_stack(columns), signal, rcond=None)[0]
    fitted = []
    for index, mode in enumerate(modes):
        cosine = coefficients[2 * index]
        sine = coefficients[2 * index + 1]
        fitted.append(
            DampedMode(
                mode.frequency_hz,
                mode.decay_seconds,
                float(np.hypot(cosine, sine)),
                float(np.arctan2(-sine, cosine)),
            )
        )
    return tuple(fitted)


def _fit_coefficients(signal: np.ndarray, poles: np.ndarray) -> np.ndarray:
    samples = np.arange(signal.size, dtype=np.float64)
    vandermonde = np.power(poles[None, :], samples[:, None])
    return np.linalg.lstsq(vandermonde, signal, rcond=None)[0]


def _mdl_value(singular_values: np.ndarray, order: int, snapshots: int) -> float:
    noise = np.maximum(
        np.square(singular_values[order:]) / snapshots, np.finfo(float).tiny
    )
    log_geometric_mean = float(np.mean(np.log(noise)))
    log_arithmetic_mean = float(np.log(np.mean(noise)))
    dimensions = singular_values.size
    likelihood = (
        -snapshots * (dimensions - order) * (log_geometric_mean - log_arithmetic_mean)
    )
    penalty = 0.5 * order * (2 * dimensions - order) * np.log(snapshots)
    return float(likelihood + penalty)


def _positive_modes(
    poles: np.ndarray, coefficients: np.ndarray, sample_rate: float
) -> list[DampedMode]:
    modes: list[DampedMode] = []
    for pole, coefficient in zip(poles, coefficients):
        angle = float(np.angle(pole))
        magnitude = float(abs(pole))
        if not 0.0 < angle < np.pi or not 0.0 < magnitude < 1.0:
            continue
        frequency = angle * sample_rate / (2.0 * np.pi)
        decay = -1.0 / (sample_rate * np.log(magnitude))
        if np.isfinite(decay) and decay > 0:
            modes.append(
                DampedMode(
                    frequency,
                    decay,
                    2.0 * abs(coefficient),
                    float(np.angle(coefficient)),
                )
            )
    return modes
