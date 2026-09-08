import numpy as np
import pytest
from triggerfish_percussion.spectral_difference import spectral_difference


def test_fixed_reference_scale_and_literal_gain():
    x = np.random.default_rng(5).normal(size=24000)
    identity = spectral_difference(x, x, 24000)
    louder = spectral_difference(x, 2 * x, 24000)
    assert np.max(np.abs(identity["difference"])) == 0
    assert identity["maximum_db"] == louder["maximum_db"]
    selected = identity["reference"] > identity["floor_db"] + 10
    assert np.allclose(louder["difference"][selected], 20 * np.log10(2))
    with pytest.raises(ValueError):
        spectral_difference(x, x[:-1], 24000)


def test_candidate_only_ringing_remains_visible():
    rate = 24000
    time = np.arange(rate) / rate
    reference = np.sin(2 * np.pi * 300 * time)
    candidate = reference + 0.2 * np.sin(2 * np.pi * 3000 * time)
    result = spectral_difference(reference, candidate, rate)
    row = np.argmin(abs(result["frequency"] - 3000))
    settled = (result["time"] > 0.2) & (result["time"] < 0.8)
    assert np.min(result["difference"][row, settled]) > 30


def test_silence_is_finite_and_neutral():
    result = spectral_difference(np.zeros(24000), np.zeros(24000), 24000)
    assert np.isfinite(result["reference"]).all()
    assert not result["difference"].any()
