"""The fast observation objective must match independently measured audio."""

from types import SimpleNamespace

import numpy as np
import pytest

from triggerfish_percussion.regional_observation_basis import RegionalObservationBasis
from triggerfish_percussion.regional_spectrum_audit import RegionalSpectrumAudit


def test_crosspower_cache_matches_welch_with_phase_and_partial_bins():
    rate = 32000
    t = np.arange(rate) / rate
    columns = np.array(
        [
            np.exp(-t * 4) * np.sin(2 * np.pi * 121.3 * t),
            np.exp(-t * 3) * np.sin(2 * np.pi * 235.7 * t + 0.6),
            np.exp(-t * 2) * np.sin(2 * np.pi * 6020.1 * t + 0.1),
        ]
    )
    amplitudes = np.array([0.3, 0.6, 0.2])
    intercept = np.random.default_rng(4).normal(0, 0.001, rate)
    original = intercept + amplitudes @ columns
    basis = SimpleNamespace(amplitudes=amplitudes, bases={7: (original, columns)})
    audit = RegionalSpectrumAudit(original, rate, [(0.03, 0.2), (0.2, 0.9)])
    cache = RegionalObservationBasis(basis, audit, 7)
    for trial in ([0.05, 0.9, 0.3], [0.7, 0.02, 0.6], amplitudes):
        trial = np.array(trial)
        assert cache.validate(basis, audit, 7, trial) < 1e-6
        actual = audit.power(intercept + trial @ columns)
        assert cache.power(trial) == pytest.approx(actual, rel=1e-7, abs=1e-12)
