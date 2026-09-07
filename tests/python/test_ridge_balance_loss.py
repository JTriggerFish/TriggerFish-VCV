"""Known false-mode and overlong-decay checks for narrow-band fitting."""

import numpy as np
from triggerfish_percussion.ridge_balance_loss import RidgeBalanceLoss


def signals():
    rate = 16000
    time = np.arange(round(1.2 * rate)) / rate
    reference = 0.5 * np.sin(2 * np.pi * 55 * time) * np.exp(-time * 30)
    reference += 0.04 * np.sin(2 * np.pi * 650 * time) * np.exp(-time * 70)
    false_mode = 0.08 * np.sin(2 * np.pi * 810 * time) * np.exp(-time * 25)
    return rate, time, reference, false_mode


def test_identical_and_false_ridge():
    rate, _, reference, false_mode = signals()
    loss = RidgeBalanceLoss(reference, reference + false_mode, rate)
    assert np.linalg.norm(loss.residual(reference)) < 1e-10
    scores = [
        np.linalg.norm(loss.ridge_residual(reference + a * false_mode, range(5)))
        for a in (0, 0.25, 0.5, 1)
    ]
    assert all(a < b for a, b in zip(scores, scores[1:]))


def test_overlong_mode_cannot_hide_in_band_total():
    rate, time, reference, _ = signals()
    excess = (
        0.04
        * np.sin(2 * np.pi * 650 * time)
        * (np.exp(-time * 20) - np.exp(-time * 70))
    )
    loss = RidgeBalanceLoss(reference, reference + excess, rate)
    assert np.linalg.norm(loss.ridge_residual(reference + excess, (2, 3))) > 1


def test_bass_ridge_and_new_quiet_bin_remain_visible():
    rate, time, reference, _ = signals()
    loss = RidgeBalanceLoss(reference, reference, rate)
    false_bass = 0.08 * np.sin(2 * np.pi * 87 * time) * np.exp(-time * 20)
    assert np.linalg.norm(loss.low.ridge_residual(reference + false_bass, range(5))) > 1
    false_high = 0.08 * np.sin(2 * np.pi * 1600 * time) * np.exp(-time * 20)
    assert np.linalg.norm(loss.ridge_residual(reference + false_high, range(5))) > 1
