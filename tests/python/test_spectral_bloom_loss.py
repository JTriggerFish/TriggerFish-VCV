"""Onset-to-bloom contrast must not be hidden by matching late decay."""

import numpy as np
import pytest

from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.spectral_bloom_basis import (
    SpectralBloomBasis,
    SpectralBloomGuard,
)
from types import SimpleNamespace
from scipy.optimize import least_squares


def sounds():
    rate = 44100
    t = np.arange(6 * rate) / rate
    low = np.sin(2 * np.pi * 180 * t) * np.exp(-2 * t)
    rng = np.random.default_rng(191)
    # A fixed, narrow random-phase high band; delayed energy, not delayed audio.
    high = sum(
        np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
        for f in np.linspace(4000, 6000, 23)
    ) / np.sqrt(23)
    rise = (1 - np.exp(-t / 0.5)) ** 4
    reference = low + 0.5 * rise * np.exp(-t) * high
    immediate = low + 0.5 * np.exp(-t) * high
    return rate, reference, immediate


def test_identity_and_gain_are_not_normalized():
    rate, reference, _ = sounds()
    loss = SpectralBloomLoss(reference, rate)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(reference * 0.5)) > 4


def test_immediate_highs_fail_even_with_correct_decay_constant():
    rate, reference, immediate = sounds()
    loss = SpectralBloomLoss(reference, rate)
    report = loss.diagnostics(immediate)
    assert report["rise_rms_db"] > 10
    assert report["envelope_rms_db"] > 5


def test_low_tone_does_not_count_as_high_frequency_bloom():
    rate, reference, _ = sounds()
    loss = SpectralBloomLoss(reference, rate)
    t = np.arange(len(reference)) / rate
    power = loss.power(np.sin(2 * np.pi * 180 * t))
    assert power[loss.edges[:-1] > 3000, 2:].max() < power.max() * 1e-8


def test_invalid_audio_rejected():
    rate, reference, _ = sounds()
    loss = SpectralBloomLoss(reference, rate)
    with pytest.raises(ValueError):
        loss.residual(reference[:-1])


def test_cached_cross_power_and_derivative():
    rate, reference, immediate = sounds()
    loss = SpectralBloomLoss(reference, rate)
    columns = np.array([reference, immediate])
    weights = np.array([0.3, 0.2])
    basis = SimpleNamespace(amplitudes=weights, bases={0: (weights @ columns, columns)})
    cache = SpectralBloomBasis(basis, loss)
    assert cache.validate(basis, np.array([0.4, 0.1])) < 1e-7
    x, direction, step = np.array([0.4, 0.1]), np.array([0.7, -0.2]), 1e-5
    numeric = (
        cache.evaluate(x + step * direction)[0]
        - cache.evaluate(x - step * direction)[0]
    ) / (2 * step)
    analytic = cache.evaluate(x)[1] @ direction
    np.testing.assert_allclose(numeric, analytic, atol=1e-6, rtol=1e-5)
    guard = SpectralBloomGuard(basis, loss, tolerance_db=1)
    assert guard.values(basis.amplitudes).min() > 0
    assert guard.values(np.array([4.0, 1.0])).min() < 0
    numeric = (
        guard.values(x + step * direction) - guard.values(x - step * direction)
    ) / (2 * step)
    np.testing.assert_allclose(
        numeric, guard.evaluate(x)[1] @ direction, atol=1e-6, rtol=1e-5
    )


def test_recovers_known_rise_and_t60_before_fitting_an_instrument():
    rate = 44100
    t = np.arange(6 * rate) / rate
    carrier = np.random.default_rng(831).normal(size=len(t))

    def render(parameters):
        rise, t60 = parameters
        return carrier * (1 - np.exp(-t / rise)) ** 4 * np.exp(-np.log(1000) * t / t60)

    loss = SpectralBloomLoss(render((0.3, 4.5)), rate)
    fitted = least_squares(
        lambda p: loss.residual(render(p)),
        (0.1, 2),
        bounds=((0.02, 0.5), (1, 12)),
        max_nfev=30,
    )
    np.testing.assert_allclose(fitted.x, (0.3, 4.5), rtol=0.01)
