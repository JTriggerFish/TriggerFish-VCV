"""Challenge the upper-layer objective with equal-power but distinct sounds."""

import numpy as np
import pytest

from triggerfish_percussion.upper_sizzle_loss import UpperSizzleLoss


@pytest.fixture(scope="module")
def signals():
    rate = 32000
    time = np.arange(6 * rate) / rate
    envelope = np.minimum(time / 0.02, 1) * np.exp(-time / 2)
    tone = lambda hz: 0.1 * envelope * np.sin(2 * np.pi * hz * time)
    return rate, time, tone


def test_upper_pitch_and_texture_are_not_hidden_by_equal_power(signals):
    rate, time, tone = signals
    plain = tone(10000)
    loss = UpperSizzleLoss(plain, rate)
    modulated = plain * (1 + 0.8 * np.sin(2 * np.pi * 40 * time)) / np.sqrt(1.32)
    assert loss.score(plain) < 1e-10
    assert loss.diagnostics(modulated)["modulation_db"] > 2
    assert loss.diagnostics(tone(13000))["centroid_octaves"] > 0.3
    assert loss.diagnostics(tone(8000))["spectrum_db"] > 3


def test_body_guard_measures_absolute_changes_without_audio_normalization(signals):
    rate, _, tone = signals
    low, high = tone(200), tone(10000)
    loss = UpperSizzleLoss(low + high, rate)
    baseline = loss.measure(low + high)
    changed = loss.diagnostics(2 * low + high, baseline)
    assert changed["low_max_change_db"] == pytest.approx(6.0206, abs=0.03)
    assert changed["low_change_db"] > 3
    assert np.asarray(baseline["low_db"]).shape == (7, 2)


def test_invalid_analysis_inputs_fail_clearly(signals):
    rate, _, tone = signals
    with pytest.raises(ValueError, match="six seconds"):
        UpperSizzleLoss(tone(10000)[:100], rate)
    loss = UpperSizzleLoss(tone(10000), rate)
    with pytest.raises(ValueError, match="finite"):
        loss.measure(np.full(6 * rate, np.nan))
