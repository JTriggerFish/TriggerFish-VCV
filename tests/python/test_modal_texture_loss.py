import numpy as np
import pytest
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss


def test_texture_distinguishes_beat_rates_without_exact_carrier_matching():
    rate = 12000
    t = np.arange(rate * 6) / rate

    def beat(frequency, spacing):
        return (
            np.sin(2 * np.pi * frequency * t)
            + np.sin(2 * np.pi * (frequency + spacing) * t)
        ) * np.exp(-t / 3)

    reference = beat(3000, 6)
    loss = ModalTextureLoss(reference, rate, centres=[3000])
    assert loss.score(reference) < 1e-12
    assert loss.score(reference * 0.25) < 1e-10
    assert loss.score(beat(3010, 6)) < loss.score(beat(3000, 35)) * 0.25
    noise = np.random.default_rng(17).normal(size=len(t)) * np.exp(-t / 3)
    assert loss.score(noise) > loss.score(beat(3010, 6)) * 4
    with pytest.raises(ValueError):
        loss.score(np.full_like(reference, np.nan))
    with pytest.raises(ValueError):
        ModalTextureLoss(np.zeros_like(reference), rate)
    with pytest.raises(ValueError):
        ModalTextureLoss(reference, rate, centres=[rate])
