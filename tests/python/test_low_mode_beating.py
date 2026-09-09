"""Known close-frequency pairs must read as slow beats, not rapid flutter."""

import numpy as np
import pytest
from triggerfish_percussion.low_mode_beating import LowModeBeating


def pair(split=1.5, t60=8):
    time = np.arange(6 * 8000) / 8000
    return (
        np.sin(2 * np.pi * 130 * time) + 0.7 * np.sin(2 * np.pi * (130 + split) * time)
    ) * 10 ** (-3 * time / t60)


def test_recovers_slow_and_fast_beating_independently_of_decay():
    loss = LowModeBeating(pair(), 8000)
    slow, fast = loss.target[0], loss.analyze(pair(12))[0]
    assert abs(slow["dominant_hz"] - 1.5) <= 0.25
    assert abs(fast["dominant_hz"] - 12) <= 0.25
    assert slow["fast_fraction"] < 0.1
    assert fast["fast_fraction"] > 0.95
    changed = loss.analyze(pair(t60=4))[0]
    np.testing.assert_allclose(changed["power"], slow["power"], rtol=0.08, atol=1e-4)
    assert loss.score(pair()) == 0


def test_single_damped_tone_does_not_look_like_beating():
    time = np.arange(6 * 8000) / 8000
    tone = np.sin(2 * np.pi * 130 * time) * 10 ** (-3 * time / 8)
    row = LowModeBeating(tone, 8000).target[0]
    assert sum(row["power"]) < 1e-5


def test_rejects_invalid_reference_and_candidate():
    with pytest.raises(ValueError):
        LowModeBeating(np.zeros(48000), 8000)
    with pytest.raises(ValueError):
        LowModeBeating(pair(), float("nan"))
    with pytest.raises(ValueError):
        LowModeBeating(pair(), 8000, region=(1, 3))
    loss = LowModeBeating(pair(), 8000)
    with pytest.raises(ValueError):
        loss.score(np.zeros(10))


def test_decay_removal_retains_beating_at_host_sample_rate():
    time = np.arange(6 * 48000) / 48000
    audio = (
        np.sin(2 * np.pi * 130 * time) + 0.7 * np.sin(2 * np.pi * 131.5 * time)
    ) * 10 ** (-3 * time / 8)
    row = LowModeBeating(audio, 48000).target[0]
    assert abs(row["dominant_hz"] - 1.5) <= 0.25
    assert row["fast_fraction"] < 0.1
