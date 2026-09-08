"""The fitting guard uses canonical measurements and exact amplitude gradients."""

from types import SimpleNamespace
import numpy as np
import pytest

torch = pytest.importorskip("torch")
from triggerfish_percussion.perceptual_envelope_guard import PerceptualEnvelopeGuard


def test_guard_measurement_and_gradient():
    torch.set_num_threads(1)
    rate = 24000
    time = np.arange(rate * 4) / rate
    reference = (
        np.random.default_rng(66).normal(size=len(time)) * 0.1 * np.exp(-4 * time)
    )
    columns = np.array(
        [reference * 0.5, np.sin(2 * np.pi * 420 * time) * np.exp(-3 * time) * 0.1]
    )
    baseline = reference * 0.8 + columns[1] * 0.1
    basis = SimpleNamespace(amplitudes=np.ones(2), bases={17: (baseline, columns)})
    guard = PerceptualEnvelopeGuard(basis, reference, rate, {17: baseline})
    assert float(guard.errors(torch.tensor(reference)).max()) < 1e-18
    assert guard.values(basis.amplitudes).min() > 0
    assert guard.values(np.array([4.0, 1.0])).min() < 0
    point = np.array([0.9, 1.03])
    _, jac = guard.evaluate(point)
    direction = np.array([0.3, -0.2])
    step = 1e-5
    numeric = (
        guard.values(point + step * direction) - guard.values(point - step * direction)
    ) / (2 * step)
    assert np.allclose(jac @ direction, numeric, rtol=1e-4, atol=1e-6)
