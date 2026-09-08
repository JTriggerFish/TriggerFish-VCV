"""Native published mel objective and its amplitude derivative agree."""

import numpy as np
import pytest
from types import SimpleNamespace

torch = pytest.importorskip("torch")
pytest.importorskip("auraloss")
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.perceptual_observation_loss import TorchAuralossMel


def test_library_score_and_amplitude_gradient():
    torch.set_num_threads(1)
    rate = 24000
    time = np.arange(rate) / rate
    reference = np.sin(2 * np.pi * 410 * time) * np.exp(-4 * time)
    reference += 0.05 * np.random.default_rng(17).normal(size=rate) * np.exp(-8 * time)
    source = AuralossMel(reference, rate)
    measurement = TorchAuralossMel(source)
    measurement.validate(reference * 0.8)
    amplitude = torch.tensor(0.8, dtype=torch.float64, requires_grad=True)
    signal = torch.tensor(reference)
    loss = measurement(signal * amplitude)
    loss.backward()
    step = 1e-5
    numeric = float(
        (measurement(signal * (0.8 + step)) - measurement(signal * (0.8 - step)))
        / (2 * step)
    )
    assert np.isclose(float(amplitude.grad), numeric, rtol=1e-4, atol=1e-6)
    # The adapter must not change the independent library scorer's precision.
    assert source.target.dtype == torch.float32
    assert source.score(reference) == 0
    from triggerfish_percussion.torch_observation_fit import _comparison_score

    search = SimpleNamespace(
        loss=source, residual=lambda _: source.residual(reference * 0.8)
    )
    assert np.isclose(
        _comparison_score(search, {}), source.score(reference * 0.8), rtol=1e-12
    )


@pytest.mark.parametrize("bounds", [(-72, 6), (-45, 7), (-40, -40), (np.nan, 6)])
def test_observation_bounds_preserve_active_ui_modes(bounds):
    from triggerfish_percussion.torch_observation_fit import polish_observation_autograd

    with pytest.raises(ValueError, match="Observation bounds"):
        polish_observation_autograd(None, bounds_db=bounds)


def test_measured_bars_can_be_frozen_without_enabling_inactive_bars():
    from triggerfish_percussion.torch_observation_fit import _observation_keys

    parameters = dict(
        resolved_level_0=-12, resolved_level_1=-6, resolved_level_2=-72, bloom_rate=1
    )
    assert _observation_keys(parameters, ("resolved_level_0",)) == ["resolved_level_1"]
    assert _observation_keys(parameters, ("resolved_level_0", "resolved_level_1")) == []
    for invalid in [("bloom_rate",), ("resolved_level_9",)]:
        with pytest.raises(ValueError):
            _observation_keys(parameters, invalid)
