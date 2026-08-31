"""Named numerical distances; no single score is an acceptance criterion."""

from dataclasses import dataclass

import numpy as np

from .modes import ModeMatching


@dataclass(frozen=True)
class LossTerm:
    name: str
    value: float
    unit: str


def energy_matching_gain(reference: np.ndarray, candidate: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=np.float64)
    current = np.asarray(candidate, dtype=np.float64)
    if (
        ref.shape != current.shape
        or not np.isfinite(ref).all()
        or not np.isfinite(current).all()
    ):
        raise ValueError("gain alignment requires matching finite arrays")
    denominator = float(np.sum(np.square(current)))
    return (
        np.sqrt(float(np.sum(np.square(ref))) / denominator) if denominator > 0 else 0.0
    )


def log_spectral_distance(
    reference_magnitude: np.ndarray,
    candidate_magnitude: np.ndarray,
    floor_db: float = -100.0,
    mask: np.ndarray | None = None,
) -> LossTerm:
    reference, candidate = _matching_arrays(reference_magnitude, candidate_magnitude)
    peak = max(float(np.max(reference)), np.finfo(float).tiny)
    floor = peak * 10.0 ** (floor_db / 20.0)
    error = 20.0 * (
        np.log10(np.maximum(candidate, floor)) - np.log10(np.maximum(reference, floor))
    )
    selected = (
        np.ones(reference.shape, dtype=bool)
        if mask is None
        else np.broadcast_to(np.asarray(mask, dtype=bool), reference.shape)
    )
    return LossTerm(
        "log_spectral", float(np.sqrt(np.mean(np.square(error[selected])))), "dB RMS"
    )


def erb_trajectory_distance(
    reference_power: np.ndarray,
    candidate_power: np.ndarray,
    floor_db: float = -100.0,
) -> tuple[LossTerm, LossTerm]:
    reference, candidate = _matching_arrays(reference_power, candidate_power)
    peak = max(float(np.max(reference)), np.finfo(float).tiny)
    floor = peak * 10.0 ** (floor_db / 10.0)
    reference_db = 10.0 * np.log10(np.maximum(reference, floor) / peak)
    candidate_db = 10.0 * np.log10(np.maximum(candidate, floor) / peak)
    level = float(np.sqrt(np.mean(np.square(candidate_db - reference_db))))
    if reference.shape[-1] > 1:
        ref_change = np.diff(reference_db, axis=-1)
        candidate_change = np.diff(candidate_db, axis=-1)
        change = float(np.sqrt(np.mean(np.square(candidate_change - ref_change))))
    else:
        change = 0.0
    return LossTerm("erb_level_trajectory", level, "dB RMS"), LossTerm(
        "erb_change_trajectory", change, "dB/frame RMS"
    )


def decay_time_distance(
    reference_seconds: np.ndarray, candidate_seconds: np.ndarray
) -> LossTerm:
    reference, candidate = _matching_arrays(reference_seconds, candidate_seconds)
    selected = (
        (reference > 0)
        & (candidate > 0)
        & np.isfinite(reference)
        & np.isfinite(candidate)
    )
    if not np.any(selected):
        return LossTerm("band_decay", float("nan"), "log-seconds RMS")
    error = np.log(candidate[selected] / reference[selected])
    return LossTerm(
        "band_decay", float(np.sqrt(np.mean(np.square(error)))), "log-seconds RMS"
    )


def modal_distance(matching: ModeMatching) -> tuple[LossTerm, ...]:
    if not matching.matches:
        nan = float("nan")
        terms = [
            LossTerm("mode_frequency", nan, "cents RMS"),
            LossTerm("mode_decay", nan, "log-seconds RMS"),
            LossTerm("mode_level", nan, "dB RMS"),
        ]
    else:
        frequency = np.array(
            [match.frequency_error_cents for match in matching.matches]
        )
        decay = np.array([match.decay_error_log for match in matching.matches])
        level = np.array([match.amplitude_error_db for match in matching.matches])
        terms = [
            LossTerm(
                "mode_frequency", float(np.sqrt(np.mean(frequency**2))), "cents RMS"
            ),
            LossTerm(
                "mode_decay", float(np.sqrt(np.mean(decay**2))), "log-seconds RMS"
            ),
            LossTerm("mode_level", float(np.sqrt(np.mean(level**2))), "dB RMS"),
        ]
    terms.append(
        LossTerm(
            "mode_count",
            float(len(matching.missing_reference) + len(matching.extra_candidate)),
            "unmatched modes",
        )
    )
    return tuple(terms)


def log_ratio_relationship_distance(
    reference_first: np.ndarray,
    reference_second: np.ndarray,
    candidate_first: np.ndarray,
    candidate_second: np.ndarray,
    name: str,
) -> LossTerm:
    ref_first, current_first = _matching_arrays(reference_first, candidate_first)
    ref_second, current_second = _matching_arrays(reference_second, candidate_second)
    if ref_first.shape != ref_second.shape:
        raise ValueError("relationship cells must share one descriptor shape")
    floor = np.finfo(float).tiny
    reference_change = np.log(
        np.maximum(ref_second, floor) / np.maximum(ref_first, floor)
    )
    candidate_change = np.log(
        np.maximum(current_second, floor) / np.maximum(current_first, floor)
    )
    value = float(np.sqrt(np.mean(np.square(candidate_change - reference_change))))
    return LossTerm(name, value, "log-ratio RMS")


def _matching_arrays(
    reference: np.ndarray, candidate: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    ref = np.asarray(reference, dtype=np.float64)
    current = np.asarray(candidate, dtype=np.float64)
    if ref.shape != current.shape or ref.size == 0:
        raise ValueError("distance inputs must have the same non-empty shape")
    if (
        np.any(ref < 0)
        or np.any(current < 0)
        or not np.isfinite(ref).all()
        or not np.isfinite(current).all()
    ):
        raise ValueError("distance inputs must be finite and non-negative")
    return ref, current
