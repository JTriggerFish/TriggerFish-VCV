"""Reject pitch, damping and gain errors before trying real drum references."""

import numpy as np
from scipy.optimize import least_squares

from triggerfish_percussion.short_drum_fit_loss import ShortDrumLoss

RATE = 16000


def hit(pitch=48, t60=0.6, phase=0):
    time = np.arange(round(1.2 * RATE)) / RATE
    return 0.3 * np.sin(2 * np.pi * pitch * time + phase) * 10 ** (-3 * time / t60)


def test_identity_gain_pitch_and_decay():
    reference = hit()
    loss = ShortDrumLoss(reference, RATE)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(reference * 2)) > 3
    assert loss.diagnostics(hit(70))["regions"][1]["low_spectrum_rmse_db"] > 5
    assert loss.diagnostics(hit(t60=0.3))["regions"][2]["envelope_rmse_db"] > 10


def test_envelope_does_not_follow_bass_cycle_phase():
    loss = ShortDrumLoss(hit(), RATE)
    # Away from the onset boundary, a quadrature sine has the same decay.
    rows = loss.diagnostics(hit(phase=np.pi / 2))["regions"]
    assert rows[2]["envelope_rmse_db"] < 0.1


def test_recover_known_pitch_and_t60_together():
    loss = ShortDrumLoss(hit(), RATE)
    result = least_squares(
        lambda x: loss.residual(hit(x[0], x[1])),
        [43, 0.8],
        bounds=([30, 0.2], [80, 1.2]),
        diff_step=0.005,
        max_nfev=30,
    )
    assert abs(result.x[0] - 48) < 0.1
    assert abs(result.x[1] - 0.6) < 0.005


def test_subsonic_energy_cannot_hide_below_the_display_range():
    reference = hit()
    loss = ShortDrumLoss(reference, RATE)
    wrong = reference + hit(10)
    assert loss.diagnostics(wrong)["regions"][1]["low_spectrum_rmse_db"] > 3
