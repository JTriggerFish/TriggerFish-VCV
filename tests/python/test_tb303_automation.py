"""Exercise the module's knob control through the production ladder/VCA DSP."""

import numpy as np
import pytest
from scipy.signal import butter, sosfilt

import _triggerfish_dsp as dsp


@pytest.mark.parametrize("sample_rate", [44100, 48000, 96000])
@pytest.mark.parametrize("update_samples", [64, 800, 1024])
def test_stepped_resonance_automation_reduces_zipper_noise(sample_rate, update_samples):
    size = sample_rate * 2
    time = np.arange(size) / sample_rate
    zero = np.zeros(size)
    # Match the report's saw oscillator, 500 Hz cutoff, and stock filter range.
    audio = dsp.tb303_oscillator_x2(zero, zero, zero, zero, zero, sample_rate)[:, 2]
    cutoff = np.full(size, 500.0)
    sweep = 0.5 + 0.45 * np.sin(2 * np.pi * 0.75 * time)
    raw = sweep[(np.arange(size) // update_samples) * update_samples]
    smoothed = dsp.tb303_resonance_control(raw, zero, sample_rate=sample_rate)
    highpass = butter(4, 8000, "highpass", fs=sample_rate, output="sos")
    noise = []
    for control in (raw, smoothed):
        rendered = dsp.diode_ladder_diagnostics_x4(
            audio, cutoff, zero, control, sample_rate=sample_rate
        )
        assert np.isfinite(rendered).all()
        assert rendered[-1, 2] == 0  # No nonlinear solver failures.
        high = sosfilt(highpass, rendered[:, 0])[sample_rate:]
        noise.append(np.sqrt(np.mean(high**2)))
    # Over 10 dB less high-frequency energy than the unsmoothed regression.
    assert noise[1] < noise[0] * 0.316


def test_resonance_cv_keeps_audio_rate_response_and_static_knob_value():
    sample_rate = 48000
    time = np.arange(sample_rate) / sample_rate
    knob = np.full(sample_rate, 0.5)
    cv = 3.0 * np.sin(2 * np.pi * 6000 * time)
    actual = dsp.tb303_resonance_control(knob, cv, amount=-0.75)
    expected = knob - 0.75 * cv / 10.0
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-15)


def test_vca_controls_affect_vca_output_and_leave_filter_output_independent():
    sample_rate = 48000
    time = np.arange(sample_rate) / sample_rate
    audio = np.sin(2 * np.pi * 220 * time)
    zero = np.zeros(sample_rate)
    cutoff = np.full(sample_rate, 500.0)
    rendered = []
    for decay in (0.0, 1.0):
        envelope = dsp.tb303_articulation(
            np.full(sample_rate, 10.0), zero, zero, vca_decay=decay
        )[:, 2]
        rendered.append(dsp.diode_ladder_voice_x4(audio, cutoff, envelope, zero))
    muted = dsp.diode_ladder_voice_x4(audio, cutoff, zero, zero)
    external = dsp.diode_ladder_voice_x4(audio, cutoff, np.ones(sample_rate), zero)
    for output in (rendered[1], muted, external):
        np.testing.assert_array_equal(output[:, 0], rendered[0][:, 0])
    assert np.max(np.abs(rendered[0][sample_rate // 2 :, 1])) < 1e-4
    assert np.sqrt(np.mean(rendered[1][sample_rate // 2 :, 1] ** 2)) > 0.01
    assert np.max(np.abs(muted[:, 1])) == 0
    assert np.sqrt(np.mean(external[:, 1] ** 2)) > 0.01
