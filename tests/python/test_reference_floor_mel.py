"""Do not reward matching quiet late hiss over a loud early spectral mismatch."""

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("auraloss")
from scipy.signal import butter, sosfiltfilt

from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.perceptual_observation_loss import TorchAuralossMel


def signal():
    rate = 24000
    t = np.arange(4 * rate) / rate
    noise = sosfiltfilt(
        butter(4, 6000, fs=rate, btype="highpass", output="sos"),
        np.random.default_rng(11).normal(size=len(t)),
    )
    high = noise * (0.12 * np.exp(-t / 0.18) + 0.00003)
    return rate, t, 0.25 * np.sin(2 * np.pi * 170 * t) * np.exp(-t / 0.8) + high, high


def test_quiet_tail_mismatch_costs_less_than_loud_early_highs():
    torch.set_num_threads(1)
    rate, t, reference, high = signal()
    loss = ReferenceFloorMel(reference, rate)
    quiet_tail_removed = reference - 0.9 * high * np.clip((t - 2) / 0.1, 0, 1)
    early_dulled = reference - 0.5 * high * np.clip((1.1 - t) / 0.1, 0, 1)
    assert loss.score(quiet_tail_removed) < 0.1 * loss.score(early_dulled)
    assert loss.score(reference) == 0
    assert loss.score(0.5 * reference) > 0.1
    floors = loss.specification["eps"].copy()
    loss.score(reference * 2)
    assert loss.specification["eps"] == floors
    # Comparison thresholds follow reference scale, even far below the
    # library's ordinary absolute STFT epsilon.
    quieter = ReferenceFloorMel(reference * 1e-6, rate)
    np.testing.assert_allclose(
        quieter.score(reference * 0.5e-6), loss.score(reference * 0.5), rtol=1e-5
    )


def test_exact_observation_autograd_adapter_retains_floor():
    torch.set_num_threads(1)
    rate, _, reference, _ = signal()
    loss = ReferenceFloorMel(reference, rate)
    adapter = TorchAuralossMel(loss)
    adapter.validate(reference * 0.8)
    gain = torch.tensor(0.8, dtype=torch.float64, requires_grad=True)
    value = adapter(torch.tensor(reference) * gain)
    value.backward()
    step = 1e-4
    derivative = (
        float(adapter(torch.tensor(reference) * (0.8 + step)))
        - float(adapter(torch.tensor(reference) * (0.8 - step)))
    ) / (2 * step)
    np.testing.assert_allclose(gain.grad.item(), derivative, rtol=1e-3, atol=1e-4)


@pytest.mark.parametrize("reference", [np.zeros(9000), np.full(9000, np.nan)])
def test_reject_invalid_reference(reference):
    with pytest.raises(ValueError):
        ReferenceFloorMel(reference, 24000)
