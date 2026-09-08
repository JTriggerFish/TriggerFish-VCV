import numpy as np
import pytest
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss


def test_unusable_reference_is_rejected():
    with pytest.raises(ValueError, match="uncontaminated"):
        BandDecayShapeLoss(np.zeros(48000 * 4), 48000)
    with pytest.raises(ValueError, match="sample rate"):
        BandDecayShapeLoss(np.zeros(48000 * 4), 0)


def test_decay_shape_is_gain_independent_but_detects_wrong_damping():
    rate = 48000
    time = np.arange(rate * 4) / rate
    noise = np.random.default_rng(11).normal(size=time.size)
    reference = noise * np.exp(-6.907755 * time / 2)
    loss = BandDecayShapeLoss(reference, rate)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(reference * 2)) < 1e-10
    wrong = noise * np.exp(-6.907755 * time / 4)
    assert np.linalg.norm(loss.residual(wrong)) > 5


def test_two_endpoint_coloured_noise_damping_is_recovered():
    from scipy.optimize import least_squares
    from scipy.signal import butter, sosfilt

    rate = 48000
    time = np.arange(rate * 5) / rate
    noise = np.random.default_rng(19).normal(size=time.size)
    bands = [(300, 700), (1500, 3000), (6000, 14000)]
    signals = [
        sosfilt(butter(5, band, btype="bandpass", fs=rate, output="sos"), noise)
        for band in bands
    ]
    positions = np.array([0, 0.5, 1])

    def render(endpoints):
        decays = np.exp(
            np.log(endpoints[0]) * (1 - positions) + np.log(endpoints[1]) * positions
        )
        return sum(
            level * signal * np.exp(-6.907755 * time / decay)
            for level, signal, decay in zip([1, 0.6, 0.3], signals, decays)
        )

    loss = BandDecayShapeLoss(render([3.5, 1.2]), rate)
    fitted = least_squares(
        lambda x: loss.residual(render(x)), [2, 2], bounds=([0.5, 0.5], [6, 6])
    )
    assert np.allclose(fitted.x, [3.5, 1.2], rtol=0.01)


def test_explicit_high_bands_are_distinct_and_gain_invariant():
    rate = 48000
    time = np.arange(rate * 4) / rate
    noise = np.random.default_rng(52).normal(size=time.size) * np.exp(-3 * time)
    bands = [(6000, 8500), (8500, 11000), (11000, 15000)]
    loss = BandDecayShapeLoss(noise, rate, bands=bands)
    assert loss.bands == tuple(bands)
    assert loss.target.shape[0] == 3
    assert np.linalg.norm(loss.residual(noise * 2)) < 1e-10
    for invalid in ([], [(100, 50)], [(0, 100)], [(100, np.nan)]):
        with pytest.raises(ValueError, match="frequency bands"):
            BandDecayShapeLoss(noise, rate, bands=invalid)
