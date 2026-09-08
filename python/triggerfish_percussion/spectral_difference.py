"""Signed spectral-power differences on a fixed reference-derived scale."""

import numpy as np
from .transforms import StftConfig, stft


def spectral_difference(reference, candidate, rate, window=8192, hop=1024):
    """Compare aligned mono audio without gain matching or time warping.

    The reference fixes the shared power floor at -70 dB. Difference pixels
    quieter than -60 dB in BOTH signals are suppressed; candidate-only energy
    remains visible. Random phase can create pixel differences even when two
    stochastic sounds share a spectrum: this is an inspection, not a loss.
    """
    if np.shape(reference) != np.shape(candidate):
        raise ValueError("Difference audio must share its duration and alignment")
    ref, model = [
        stft(x, rate, StftConfig(window, hop)) for x in (reference, candidate)
    ]
    maximum = max(float(ref.power.max()), 1e-20)
    floor = maximum * 1e-7
    ref_db = 10 * np.log10(np.maximum(ref.power, floor))
    model_db = 10 * np.log10(np.maximum(model.power, floor))
    visible = (ref.power > floor * 10) | (model.power > floor * 10)
    return dict(
        frequency=ref.frequencies_hz,
        time=ref.times_seconds,
        reference=ref_db,
        candidate=model_db,
        difference=np.where(visible, model_db - ref_db, 0),
        maximum_db=10 * np.log10(maximum),
        floor_db=10 * np.log10(floor),
    )
