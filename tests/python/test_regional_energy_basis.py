from types import SimpleNamespace
import numpy as np
from triggerfish_percussion.regional_energy_basis import (
    RegionalEnergyBasis,
    RegionalEnergyGuard,
)
from triggerfish_percussion.regional_energy_loss import RegionalEnergyLoss


def test_correlated_power_and_gradient():
    rate = 24000
    time = np.arange(rate) / rate
    first = np.sin(2 * np.pi * 410 * time) * np.exp(-3 * time)
    second = 0.8 * first + 0.05 * np.random.default_rng(17).normal(size=rate)
    columns = np.array([first, second])
    amplitudes = np.array([0.3, 0.7])
    baseline = amplitudes @ columns + 0.03 * np.cos(2 * np.pi * 220 * time)
    basis = SimpleNamespace(
        amplitudes=amplitudes,
        bases={1: (baseline, columns), 2: (baseline * 0.9, columns)},
    )
    loss = RegionalEnergyLoss(
        baseline, rate, [(100, 1000), (1000, 8000)], [(0, 0.1), (0.1, 0.5), (0.5, 1)]
    )
    prepared = RegionalEnergyBasis(basis, loss)
    point = np.array([0.9, 0.2])
    assert prepared.validate(basis, point) < 1e-10
    direction = np.array([0.3, -0.1])
    step = 1e-5
    numeric = (
        prepared.evaluate(point + step * direction)[0]
        - prepared.evaluate(point - step * direction)[0]
    ) / (2 * step)
    assert np.allclose(
        prepared.evaluate(point)[1] @ direction, numeric, rtol=1e-5, atol=1e-7
    )
    guard = RegionalEnergyGuard(basis, loss, {1: baseline, 2: baseline * 0.9})
    assert guard.values(amplitudes).min() > 0
    assert guard.values(np.array([4.0, 4.0])).min() < 0
    numeric = (
        guard.values(point + step * direction) - guard.values(point - step * direction)
    ) / (2 * step)
    assert np.allclose(
        guard.evaluate(point)[1] @ direction, numeric, rtol=1e-5, atol=1e-7
    )
