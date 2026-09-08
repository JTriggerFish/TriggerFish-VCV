import numpy as np
import pytest
from triggerfish_percussion.regional_energy_loss import RegionalEnergyLoss


def fixture():
    rate = 24000
    audio = np.random.default_rng(12).normal(size=rate) * 0.1
    loss = RegionalEnergyLoss(
        audio, rate, [(100, 1000), (1000, 8000)], [(0, 0.1), (0.1, 0.5), (0.5, 1)]
    )
    return audio, loss


def test_identity_and_absolute_gain():
    audio, loss = fixture()
    assert np.max(np.abs(loss.residual(audio))) == 0
    assert np.allclose(loss.diagnostics(audio * 2)["difference_db"], 20 * np.log10(2))


def test_late_gain_does_not_renormalize_early_energy():
    audio, loss = fixture()
    changed = audio.copy()
    changed[12000:] *= 0.5
    error = np.array(loss.diagnostics(changed)["difference_db"])
    assert np.all(error[:, :2] == 0)
    assert np.all(error[:, 2] < -5.9)


def test_invalid_audio_and_regions():
    audio, loss = fixture()
    with pytest.raises(ValueError):
        loss.residual(audio[:-1])
    with pytest.raises(ValueError):
        RegionalEnergyLoss(audio, 24000, [(100, 1000)], [(0, 0.6), (0.5, 1)])
