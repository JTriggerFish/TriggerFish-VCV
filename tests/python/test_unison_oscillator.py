import numpy as np

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000


def render(*, voices=7, spread=4.0, pulse_mix=0.0, width=0.65, seconds=1.0):
    samples = round(SAMPLE_RATE * seconds)
    frequency = np.full(samples, 440.0)
    pulse_width = 0.5 + 0.4 * np.sin(
        2.0 * np.pi * 0.7 * np.arange(samples) / SAMPLE_RATE
    )
    return dsp.stacked_oscillator(
        frequency,
        pulse_width,
        voices=voices,
        spread_cents=spread,
        pulse_mix=pulse_mix,
        width=width,
        sample_rate=SAMPLE_RATE,
    )


def test_large_stack_layouts_are_centered_and_bounded():
    for voices in range(1, 17):
        pitch = dsp.stacked_oscillator_pitch_positions(voices)
        pan = dsp.stacked_oscillator_pan_positions(voices)
        assert abs(np.mean(pitch)) < 1e-14
        assert abs(np.mean(pan)) < 1e-14
        assert np.max(np.abs(pitch)) <= 1.0
        assert np.max(np.abs(pan)) <= 1.0


def test_stacked_saw_and_pwm_pulse_are_finite():
    for pulse_mix in (0.0, 0.35, 1.0):
        output = render(voices=16, pulse_mix=pulse_mix, seconds=0.25)
        assert output.shape == (SAMPLE_RATE // 4, 4)
        assert np.all(np.isfinite(output))
        assert np.max(np.abs(output)) < 20.0


def test_zero_width_stereo_matches_mono_contract():
    output = render(voices=9, width=0.0, seconds=0.25)
    np.testing.assert_allclose(output[:, 1], output[:, 0], atol=2e-12)
    np.testing.assert_allclose(output[:, 2], output[:, 0], atol=2e-12)


def test_default_detuned_stack_has_rough_energy_normalization():
    warmup = SAMPLE_RATE // 2
    single = render(voices=1, seconds=2.0)[warmup:, 0]
    for voices in (4, 7, 16):
        stacked = render(voices=voices, seconds=2.0)[warmup:, 0]
        ratio = np.sqrt(np.mean(stacked**2) / np.mean(single**2))
        assert 0.75 < ratio < 1.35
