"""Dedicated T60 acceptance for crash-cymbal calibration."""

from __future__ import annotations

import numpy as np

from .alignment import detect_impact_onset
from .t60_envelope import fit_two_point_t60, measure_band_t60

MAXIMUM_BAND_T60_RMSE_LOG = 0.3
MAXIMUM_BAND_T60_ABSOLUTE_LOG = 0.6
MINIMUM_COMMON_T60_BANDS = 4
MAXIMUM_BODY_T60_FREQUENCY_HZ = 15_000.0
MINIMUM_T60_R_SQUARED = 0.95
T60_START_SECONDS = 0.05


def t60_diagnostics(reference, candidate) -> dict[str, object]:
    """Compare only decay bands that the reference recording can support."""
    reference_fit = _measured_fit(
        reference.samples,
        reference.sample_rate,
        detect_impact_onset(reference.samples, reference.sample_rate),
    )
    if reference_fit is None or (
        reference_fit.band_frequencies_hz.size < MINIMUM_COMMON_T60_BANDS
    ):
        return {
            "passed": False,
            "evaluated": False,
            "status": "insufficient_reference_evidence",
            "reference_band_count": (
                0
                if reference_fit is None
                else int(reference_fit.band_frequencies_hz.size)
            ),
        }
    candidate_fit = _measured_fit(candidate.samples, candidate.sample_rate, 0)
    if candidate_fit is None:
        return {
            "passed": False,
            "evaluated": True,
            "status": "insufficient_candidate_evidence",
            "reference_band_count": int(reference_fit.band_frequencies_hz.size),
            "common_band_count": 0,
        }

    endpoint_difference = np.log(
        [candidate_fit.low_seconds, candidate_fit.high_seconds]
    ) - np.log([reference_fit.low_seconds, reference_fit.high_seconds])
    endpoint_rmse = float(np.sqrt(np.mean(np.square(endpoint_difference))))
    endpoint_maximum = float(np.max(np.abs(endpoint_difference)))
    frequencies, reference_bands, candidate_bands = _common_bands(
        reference_fit, candidate_fit
    )
    band_difference = np.log(candidate_bands / reference_bands)
    band_rmse = _rmse_or_infinity(band_difference)
    band_maximum = _maximum_or_infinity(band_difference)
    enough_common_bands = band_difference.size >= MINIMUM_COMMON_T60_BANDS
    return {
        "passed": bool(
            enough_common_bands
            and band_rmse <= MAXIMUM_BAND_T60_RMSE_LOG
            and band_maximum <= MAXIMUM_BAND_T60_ABSOLUTE_LOG
        ),
        "evaluated": True,
        "status": (
            "measured" if enough_common_bands else "insufficient_candidate_evidence"
        ),
        "reference_endpoint_seconds": [
            reference_fit.low_seconds,
            reference_fit.high_seconds,
        ],
        "candidate_endpoint_seconds": [
            candidate_fit.low_seconds,
            candidate_fit.high_seconds,
        ],
        "reference_curve_log_rmse": reference_fit.log_rmse,
        "candidate_curve_log_rmse": candidate_fit.log_rmse,
        "endpoint_rmse_log": endpoint_rmse,
        "endpoint_maximum_absolute_log": endpoint_maximum,
        "common_band_count": int(band_difference.size),
        "band_rmse_log": band_rmse,
        "band_maximum_absolute_log": band_maximum,
        "band_frequencies_hz": frequencies.tolist(),
        "reference_band_seconds": reference_bands.tolist(),
        "candidate_band_seconds": candidate_bands.tolist(),
    }


def _measured_fit(samples, sample_rate, onset_sample):
    try:
        frequencies, fits = measure_band_t60(
            samples,
            sample_rate,
            onset_sample=onset_sample,
            start_seconds=T60_START_SECONDS,
            maximum_hz=MAXIMUM_BODY_T60_FREQUENCY_HZ,
            band_count=24,
        )
        return fit_two_point_t60(
            frequencies,
            fits,
            sample_rate,
            minimum_r_squared=MINIMUM_T60_R_SQUARED,
        )
    except ValueError:
        return None


def _common_bands(reference, candidate):
    reference_indices = _frequency_indices(reference.band_frequencies_hz)
    candidate_indices = _frequency_indices(candidate.band_frequencies_hz)
    common = sorted(reference_indices.keys() & candidate_indices.keys())
    return (
        np.asarray(common),
        np.asarray(
            [reference.measured_seconds[reference_indices[key]] for key in common]
        ),
        np.asarray(
            [candidate.measured_seconds[candidate_indices[key]] for key in common]
        ),
    )


def _frequency_indices(frequencies) -> dict[float, int]:
    return {
        round(float(frequency), 6): index for index, frequency in enumerate(frequencies)
    }


def _rmse_or_infinity(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values)))) if values.size else float("inf")


def _maximum_or_infinity(values: np.ndarray) -> float:
    return float(np.max(np.abs(values))) if values.size else float("inf")
