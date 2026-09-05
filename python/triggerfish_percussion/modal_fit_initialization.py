"""Spectral mode-placement proposals, not claims of physical modal decomposition."""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from .transforms import StftConfig, stft


def spectral_mode_candidates(
    samples, sample_rate, *, count=16, low=25, high=8000, start=0.03, end=0.35
):
    """Rank resolved peaks across log bands; broadband energy supplies no fake peaks.

    Long windows resolve bass but smear attack timing. These frequencies only
    seed full-render comparisons; their measured powers are not modal gains.
    """
    if count < 1 or not 0 < low < high or not 0 <= start < end:
        raise ValueError("Invalid modal proposal ranges")
    transformed = stft(samples, sample_rate, StftConfig(8192, 256))
    selected = (transformed.times_seconds >= start) & (transformed.times_seconds < end)
    if not selected.any():
        raise ValueError("No reference frames inside the modal proposal window")
    power = gaussian_filter1d(transformed.power[:, selected].mean(axis=1), 0.6)
    db = 10 * np.log10(np.maximum(power, max(float(power.max()) * 1e-10, 1e-30)))
    peaks, info = find_peaks(db, prominence=3)
    candidates = [
        dict(
            frequency=float(transformed.frequencies_hz[i]),
            power_db=float(db[i]),
            prominence_db=float(prominence),
        )
        for i, prominence in zip(peaks, info["prominences"])
        if low <= transformed.frequencies_hz[i] <= high and db[i] >= db.max() - 55
    ]
    ranked = sorted(candidates, key=lambda item: item["power_db"], reverse=True)
    # Reserve room for quieter upper ringing; do not let one bass cluster take
    # every slot. Remaining slots still go to the strongest resolved peaks.
    chosen = []
    for left, right in zip(
        np.geomspace(low, high, 5)[:-1], np.geomspace(low, high, 5)[1:]
    ):
        band = [item for item in ranked if left <= item["frequency"] <= right]
        if band and len(chosen) < count:
            chosen.append(band[0])
    for item in ranked:
        if len(chosen) >= count:
            break
        if item not in chosen:
            chosen.append(item)
    return sorted(chosen, key=lambda item: item["frequency"])


def reference_modal_starts(
    parameters, proposals, *, prefix="resonance", capacity=16, counts=(4, 8, 12, 16)
):
    """Write explicit parameter vectors; no runtime formula or new fitted controls.

    Prominence uses a few neutral falloffs instead of equating measured spectrum
    with modal coefficients (which also depend on excitation and observation).
    Spatial couplings are unit starting values, visible in each saved mode.
    """
    if not proposals:
        return []
    ranked = sorted(proposals, key=lambda item: item["power_db"], reverse=True)
    results = []
    for count in sorted({min(n, capacity, len(proposals)) for n in counts}):
        selected = sorted(ranked[:count], key=lambda item: item["frequency"])
        for falloff in (0, 6, 12):
            values = dict(parameters)
            for i in range(capacity):
                values[f"{prefix}_level_{i}"] = -72
                if i < count:
                    frequency = selected[i]["frequency"]
                    values[f"{prefix}_frequency_{i}"] = frequency
                    values[f"{prefix}_level_{i}"] = max(
                        -60, -falloff * np.log2(frequency / selected[0]["frequency"])
                    )
                    values[f"{prefix}_centre_{i}"] = values[f"{prefix}_edge_{i}"] = 1
            results.append(
                (f"reference peaks: {count} modes, {falloff} dB/oct", values)
            )
    return results
