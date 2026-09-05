import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from triggerfish_percussion.power_envelope import smoothed_power


@pytest.mark.parametrize("sigma", [0.3, 4, 192, 529.2, 1152])
def test_matches_direct_reflected_gaussian(sigma):
    samples = np.random.default_rng(2).normal(size=2048)
    samples[0], samples[-1] = 10, -10
    np.testing.assert_allclose(
        smoothed_power(samples, sigma),
        gaussian_filter1d(samples**2, sigma),
        atol=1e-12,
        rtol=1e-12,
    )


def test_silence_and_decay_floor():
    assert np.all(smoothed_power(np.zeros(1024), 200) == 0)
    samples = np.exp(-np.arange(8192) / 60)
    actual = smoothed_power(samples, 200)
    expected = gaussian_filter1d(samples**2, 200)
    np.testing.assert_allclose(actual, expected, atol=1e-14)
