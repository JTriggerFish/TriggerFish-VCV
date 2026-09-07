"""Band balance must expose weak bass and an extended noisy attack."""

import numpy as np
from triggerfish_percussion.drum_balance_loss import DrumBalanceLoss

RATE = 16000


def sources():
    time = np.arange(round(1.2 * RATE)) / RATE
    bass = 0.4 * np.sin(2 * np.pi * 65 * time) * 10 ** (-3 * time / 0.35)
    noise = np.random.default_rng(12).normal(size=len(time))
    return time, bass, noise


def test_identity_and_no_candidate_normalization():
    time, bass, noise = sources()
    reference = bass + 0.03 * noise * 10 ** (-3 * time / 0.03)
    loss = DrumBalanceLoss(reference, RATE)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.balance_residual(reference * 2)) > 2


def test_weak_bass_and_long_noise_are_penalized():
    time, bass, noise = sources()
    attack = 0.03 * noise * 10 ** (-3 * time / 0.03)
    loss = DrumBalanceLoss(bass + attack, RATE)
    assert np.linalg.norm(loss.balance_residual(0.5 * bass + attack)) > 2
    wrong = bass + 0.03 * noise * 10 ** (-3 * time / 0.2)
    assert np.linalg.norm(loss.balance_residual(wrong)) > 2


def test_silence_and_low_sample_rate_are_finite():
    silence = np.zeros(9600)
    loss = DrumBalanceLoss(silence, 8000)
    assert np.all(np.isfinite(loss.residual(silence)))
