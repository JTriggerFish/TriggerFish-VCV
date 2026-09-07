"""Regression fixtures for failure modes observed across metallic instruments."""

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar

from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss


@pytest.fixture(scope="module")
def example():
    rate = 24000
    t = np.arange(rate * 6) / rate
    rng = np.random.default_rng(93)
    noise = 0.025 * rng.normal(size=len(t)) * np.exp(-t / 0.8)
    tone = 0.15 * np.sin(2 * np.pi * 515 * t) * np.exp(-t / 1.2)
    tone += 0.08 * np.sin(2 * np.pi * 1680 * t) * np.exp(-t / 0.7)
    return t, tone, noise, MetallicBalanceLoss(tone + noise, rate)


def test_identity_and_gain_recovery(example):
    _, tone, noise, loss = example
    x = tone + noise
    assert np.linalg.norm(loss.residual(x)) == 0
    scores = [np.linalg.norm(loss.residual(x * g)) for g in (1, 0.5, 0.25)]
    assert scores[0] < scores[1] < scores[2]
    result = minimize_scalar(
        lambda gain: np.linalg.norm(loss.residual(x * gain)),
        bounds=(0.1, 2),
        method="bounded",
        options={"xatol": 0.005},
    )
    assert abs(result.x - 1) < 0.01


def test_missing_attack_and_ringing_have_separate_diagnostics(example):
    t, tone, noise, loss = example
    attack = (tone + noise).copy()
    attack[t < 0.03] *= 0.1
    parts = loss.diagnostics(attack)["components"]
    assert parts["envelope"] > 1
    assert parts["attack"] > 2
    washed = gaussian_filter1d(tone, 2) + noise * 2
    assert loss.diagnostics(washed)["components"]["contrast"] > 0.5
    assert loss.diagnostics(noise)["components"]["linear_spectrum"] > 10


def test_late_tail_is_not_ignored():
    rate = 24000
    t = np.arange(rate * 10) / rate
    reference = 0.2 * np.sin(2 * np.pi * 430 * t) * np.exp(-t / 4)
    loss = MetallicBalanceLoss(reference, rate)
    clipped = reference.copy()
    clipped[t > 6] = 0
    assert np.linalg.norm(loss.residual(clipped, (4,))) > 5
    assert loss.features.regions[-1] == (3, 10)


def test_shifted_bloom_is_detected():
    rate = 24000
    t = np.arange(rate * 6) / rate
    rng = np.random.default_rng(12)
    carrier = rng.normal(size=len(t))
    envelope = t * np.exp(-t / 0.6)
    ref = envelope * carrier
    loss = MetallicBalanceLoss(ref, rate)
    early = np.exp(-t / 0.6) * carrier
    early *= np.linalg.norm(ref) / np.linalg.norm(early)
    assert loss.diagnostics(early)["components"]["envelope"] > 5


def test_mismatched_duration_rejected(example):
    _, tone, noise, loss = example
    with pytest.raises(ValueError, match="duration"):
        loss.residual((tone + noise)[:-10])


def test_empty_resolved_region_rejected():
    # A 4096-sample hop at 8 kHz skips the entire 120–500 ms region.
    with pytest.raises(ValueError, match="no frames in region"):
        MetallicBalanceLoss(np.zeros(8000 * 4), 8000)


def test_complex_audio_rejected(example):
    _, tone, _, loss = example
    with pytest.raises(ValueError, match="real mono"):
        loss.residual(tone.astype(complex))


def test_missing_first_ten_ms_not_hidden_by_matched_total_energy(example):
    t, tone, noise, _ = example
    reference = tone + noise
    altered = reference.copy()
    altered[t < 0.01] *= 0.1
    region = (t >= 0.01) & (t < 0.1)
    lost = np.sum(reference[t < 0.01] ** 2) - np.sum(altered[t < 0.01] ** 2)
    altered[region] *= np.sqrt(1 + lost / np.sum(altered[region] ** 2))
    loss = MetallicBalanceLoss(reference, 24000, fast_attack=True)
    assert np.sum(reference[t < 0.1] ** 2) == pytest.approx(
        np.sum(altered[t < 0.1] ** 2)
    )
    assert loss.diagnostics(altered)["components"]["attack"] > 10
