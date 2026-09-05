"""Recovery checks before trusting a stochastic percussion search objective."""

import numpy as np
from scipy.optimize import minimize_scalar

from triggerfish_percussion.trajectory_fit_loss import TrajectoryLoss


def fixture(t60=3, pitch=375):
    sample_rate = 24000
    time = np.arange(6 * sample_rate) / sample_rate
    noise = np.random.default_rng(19).normal(size=len(time))
    source = 0.3 * np.sin(2 * np.pi * pitch * time) + 0.08 * noise
    return sample_rate, source * 10 ** (-3 * time / t60)


def test_identical_signal_and_literal_gain():
    # Stay above the reference floor throughout all five measured regions.
    rate, reference = fixture(t60=20)
    loss = TrajectoryLoss(reference, rate)
    assert np.max(np.abs(loss.residual(reference))) == 0
    # A 6 dB level error must not vanish through candidate normalization.
    assert 5.5 < np.linalg.norm(loss.residual(reference * 2)) < 6.1


def test_wrong_pitch_is_penalized():
    rate, reference = fixture()
    loss = TrajectoryLoss(reference, rate)
    wrong = fixture(pitch=550)[1]
    diagnostics = loss.diagnostics(wrong)
    assert diagnostics["regions"][0]["ridge_rmse_db"] > 5


def test_known_coloured_decay_is_recovered():
    rate, reference = fixture(t60=4.2)
    loss = TrajectoryLoss(reference, rate)
    result = minimize_scalar(
        lambda t60: np.linalg.norm(loss.residual(fixture(t60=t60)[1])),
        bounds=(1, 8),
        method="bounded",
        options={"xatol": 0.02},
    )
    assert abs(result.x - 4.2) < 0.03


def test_late_mismatch_cannot_hide_behind_correct_attack():
    rate, reference = fixture(t60=4.2)
    loss = TrajectoryLoss(reference, rate)
    wrong = reference.copy()
    wrong[rate:] *= 0.1
    rows = loss.diagnostics(wrong)["regions"]
    assert rows[0]["band_rmse_db"] < 0.01
    assert rows[3]["band_rmse_db"] > 15
