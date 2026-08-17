import numpy as np
import pytest

import _triggerfish_dsp as dsp

# Continuous-time limit-cycle frequencies generated with SciPy's fifth-order
# Radau solver at rtol=1e-9 and atol=1e-11.
REFERENCE_FREQUENCY_SCALE = {
    0.1: 0.9993755712,
    0.5: 0.9847209766,
    2.0: 0.8234978601,
    5.0: 0.5410834048,
    9.0: 0.3579717065,
}


def render_vdpo(frequency, damping, sample_rate=48_000, duration=0.45):
    size = int(sample_rate * duration)
    return dsp.vdpo(
        np.zeros(size),
        np.full(size, damping),
        np.full(size, 2.0 * np.pi * frequency),
        sample_rate,
    )


def render_vdpo_pitch(processor, frequency, damping, sample_rate=48_000, duration=0.8):
    size = int(sample_rate * duration)
    return processor(
        np.zeros(size),
        np.full(size, damping),
        np.full(size, np.log2(frequency / 261.625565)),
        sample_rate,
    )


def estimate_frequency(signal, sample_rate=48_000, settle_time=0.15):
    settled = signal[int(sample_rate * settle_time) :]
    indices = np.flatnonzero((settled[:-1] <= 0.0) & (settled[1:] > 0.0))
    assert indices.size >= 3
    crossings = indices - settled[indices] / (settled[indices + 1] - settled[indices])
    return (crossings.size - 1) * sample_rate / (crossings[-1] - crossings[0])


def test_vca_binding_is_finite_and_length_preserving():
    audio = np.sin(np.linspace(0.0, 8.0 * np.pi, 1024, dtype=np.float32))
    cv = np.linspace(0.0, 1.0, audio.size, dtype=np.float32)
    output = dsp.vca_transistor(audio, cv, 48_000.0)
    assert output.shape == audio.shape
    assert np.isfinite(output).all()
    assert np.max(np.abs(output)) < 11.7


def test_slop_linear_detune_is_exact_when_drift_is_zero():
    pitch = np.linspace(-10.0, 10.0, 4097)
    output = dsp.detune_optimized(pitch, np.zeros_like(pitch))
    assert np.array_equal(output, pitch)


def test_slop_linear_detune_stays_within_musical_pitch_error_budget():
    pitch = np.repeat(np.linspace(-5.0, 5.0, 401), 401)
    detune = np.tile(np.linspace(-2.0, 2.0, 401), 401)
    reference = dsp.detune_reference_double(pitch, detune)
    optimized = dsp.detune_optimized(pitch, detune)
    error_cents = 1200.0 * (optimized - reference)
    assert np.max(np.abs(error_cents)) < 0.002


def test_slop_linear_detune_preserves_wide_range_precision():
    pitch = np.repeat(np.linspace(-10.0, 10.0, 401), 101)
    relative_detune = np.tile(np.linspace(-0.25, 0.25, 101), 401)
    detune = 261.63 * np.exp2(pitch) * relative_detune
    reference = dsp.detune_reference_double(pitch, detune)
    optimized = dsp.detune_optimized(pitch, detune)
    error_cents = 1200.0 * (optimized - reference)
    assert np.max(np.abs(error_cents)) < 0.004


def test_vdpo_binding_is_finite_and_length_preserving():
    size = 256
    audio = np.zeros(size)
    damping = np.full(size, 0.5)
    angular_frequency = np.full(size, 2.0 * np.pi * 261.625565)
    output = dsp.vdpo(audio, damping, angular_frequency, 48_000.0)
    assert output.shape == audio.shape
    assert np.isfinite(output).all()
    assert np.max(np.abs(output)) < 11.7


@pytest.mark.parametrize("damping", (0.5, 9.0))
@pytest.mark.parametrize("frequency", (20.0, 1_000.0, 8_000.0, 18_000.0))
def test_shared_fast_exp2_stays_within_pitch_error_budget(damping, frequency):
    standard_frequency = estimate_frequency(
        render_vdpo_pitch(dsp.vdpo_pitch_std, frequency, damping)
    )
    fast_frequency = estimate_frequency(
        render_vdpo_pitch(dsp.vdpo_pitch_fast_exp2, frequency, damping)
    )
    pitch_error_cents = 1200.0 * np.log2(fast_frequency / standard_frequency)
    assert abs(pitch_error_cents) < 0.02


@pytest.mark.parametrize("damping, expected_scale", REFERENCE_FREQUENCY_SCALE.items())
def test_vdpo_matches_high_precision_continuous_period(damping, expected_scale):
    requested_frequency = 523.25113
    measured_frequency = estimate_frequency(render_vdpo(requested_frequency, damping))
    assert measured_frequency / requested_frequency == pytest.approx(
        expected_scale, rel=2.0e-3
    )


@pytest.mark.parametrize("damping", REFERENCE_FREQUENCY_SCALE)
def test_vdpo_tracks_pitch_above_legacy_limit(damping):
    reference_frequency = 523.25113
    reference_scale = (
        estimate_frequency(render_vdpo(reference_frequency, damping))
        / reference_frequency
    )

    for requested_frequency in (1_046.502, 4_186.009, 8_000.0, 12_000.0, 18_000.0):
        measured_frequency = estimate_frequency(
            render_vdpo(requested_frequency, damping)
        )
        tracking_scale = measured_frequency / requested_frequency
        assert tracking_scale == pytest.approx(reference_scale, rel=1.2e-2)


@pytest.mark.parametrize("sample_rate", (44_100, 96_000))
@pytest.mark.parametrize("damping", (0.5, 2.0, 9.0))
def test_vdpo_tracking_is_sample_rate_independent(sample_rate, damping):
    reference_frequency = 523.25113
    reference_scale = (
        estimate_frequency(
            render_vdpo(reference_frequency, damping, sample_rate), sample_rate
        )
        / reference_frequency
    )
    requested_frequency = 0.4 * sample_rate
    measured_frequency = estimate_frequency(
        render_vdpo(requested_frequency, damping, sample_rate), sample_rate
    )
    assert measured_frequency / requested_frequency == pytest.approx(
        reference_scale, rel=1.2e-2
    )


def test_vdpo_remains_finite_during_audio_rate_parameter_sweep():
    size = 16_384
    sample_rate = 48_000.0
    output = dsp.vdpo(
        np.linspace(-5.0, 5.0, size),
        np.linspace(0.001, 9.0, size),
        2.0 * np.pi * np.geomspace(20.0, 21_600.0, size),
        sample_rate,
    )
    assert np.isfinite(output).all()
    assert np.max(np.abs(output)) < 11.7


def test_bindings_reject_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        dsp.vca_transistor(np.zeros(4), np.zeros(3), 48_000.0)
