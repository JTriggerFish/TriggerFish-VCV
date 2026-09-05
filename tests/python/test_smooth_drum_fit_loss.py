"""Smooth objective remains level-sensitive and recovers bass pitch/decay."""

import numpy as np
from scipy.optimize import least_squares

from triggerfish_percussion.smooth_drum_fit_loss import SmoothDrumLoss


def fixture(pitch=48, t60=0.6):
    time = np.arange(19200) / 16000
    return 0.3 * np.sin(2 * np.pi * pitch * time) * 10 ** (-3 * time / t60)


def test_identity_and_no_level_matching():
    reference = fixture()
    loss = SmoothDrumLoss(reference, 16000)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(2 * reference)) > 0.05


def test_recover_bass_from_different_starting_pitches():
    loss = SmoothDrumLoss(fixture(), 16000)
    for pitch in (35, 43, 65):
        result = least_squares(
            lambda x: loss.residual(fixture(x[0], x[1])),
            [pitch, 0.9],
            bounds=([25, 0.2], [80, 1.2]),
            jac="3-point",
            diff_step=0.0001,
            max_nfev=80,
        )
        assert abs(result.x[0] - 48) < 0.15
        assert abs(result.x[1] - 0.6) < 0.01
