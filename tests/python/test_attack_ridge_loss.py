import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("auraloss")
from triggerfish_percussion.attack_ridge_loss import AttackRidgeLoss


def test_attack_ridges_distinguish_small_pitch_errors_and_fixed_gain():
    torch.set_num_threads(1)
    rate = 44100
    time = np.arange(rate) / rate

    def tone(hz):
        return np.sin(2 * np.pi * hz * time) * np.exp(-time * 3)

    reference = tone(123.47)
    loss = AttackRidgeLoss(reference, rate)
    assert loss.score(reference) < 1e-6
    assert loss.score(tone(127)) > 0.05
    assert loss.score(reference * 0.5) > 0.1
    changed_tail = reference.copy()
    changed_tail[round(0.4 * rate) :] = 0
    assert loss.score(changed_tail) < 1e-6
    with pytest.raises(ValueError):
        loss.score(reference[:100])
    with pytest.raises(ValueError):
        loss.score(np.full(rate, np.nan))
    with pytest.raises(ValueError):
        AttackRidgeLoss(np.full(rate, np.nan), rate)
    with pytest.raises(ValueError):
        AttackRidgeLoss(reference, 0)
