import numpy as np
import pytest

import _triggerfish_dsp as dsp
from reference_arp4072 import (
    circuit_values,
    harmonic_amplitudes,
    midpoint_transfer,
    render_continuous_sine,
)

SAMPLE_RATE = 48_000.0


def _sine_gain(frequency, cutoff, resonance=0.0, amplitude=1.0e-3):
    duration = max(1.0, 80.0 / frequency)
    sample_count = int(SAMPLE_RATE * duration)
    time = np.arange(sample_count) / SAMPLE_RATE
    signal = amplitude * np.sin(2.0 * np.pi * frequency * time)
    # The x1 path isolates the prewarped solver from resampler group delay.
    output = dsp.arp4072_x1(signal, cutoff, resonance)
    settle = sample_count // 2
    reference = np.exp(-2j * np.pi * frequency * time[settle:])
    phasor = 2.0 * np.mean(output[settle:] * reference)
    # A sine's positive-frequency phasor is -j times its peak amplitude.
    return phasor / (-1j * amplitude)


def test_circuit_ratios_match_the_independent_reference():
    expected = circuit_values()
    actual = dsp.arp4072_circuit_values()
    assert actual.keys() == expected.keys()
    for name, value in expected.items():
        assert np.isclose(actual[name], value, rtol=1.0e-13, atol=1.0e-15), name

    ratio = actual["feedback_base_scale"] / actual["audio_base_scale"]
    assert 6.57 < ratio < 6.60
    assert 0.98 < actual["small_signal_input_gain"] < 1.03
    assert np.isclose(actual["stage_base_resistance_ohms"], 212.280701754, rtol=1e-11)


def test_small_signal_response_matches_the_circuit_transfer():
    cutoff = 1_000.0
    resonance = 0.35
    for frequency in (50.0, 300.0, 1_000.0, 3_000.0):
        measured = _sine_gain(frequency, cutoff, resonance)
        expected = midpoint_transfer(frequency, cutoff, resonance, SAMPLE_RATE)
        magnitude_error_db = 20.0 * np.log10(abs(measured / expected))
        phase_error = np.angle(measured / expected)
        assert abs(magnitude_error_db) < 0.12
        assert abs(phase_error) < 0.025


def test_input_pair_compresses_smoothly_without_a_hard_boundary():
    frequency = 100.0
    cutoff = 8_000.0
    sample_count = int(SAMPLE_RATE)
    time = np.arange(sample_count) / SAMPLE_RATE

    gains = []
    maximum_steps = []
    for amplitude in (0.1, 1.0, 5.0, 10.0):
        signal = amplitude * np.sin(2.0 * np.pi * frequency * time)
        output = dsp.arp4072_x4(signal, cutoff, 0.0)
        steady = output[sample_count // 2 :]
        gains.append(np.sqrt(np.mean(steady**2)) / (amplitude / np.sqrt(2.0)))
        maximum_steps.append(np.max(np.abs(np.diff(steady))))
        assert np.all(np.isfinite(output))

    assert gains[0] > gains[1] > gains[2] > gains[3]
    assert maximum_steps[-1] < 0.25


def test_high_resonance_self_oscillation_is_bounded_and_smooth():
    silence = np.zeros(int(2.0 * SAMPLE_RATE))
    silence[0] = 1.0e-6
    output = dsp.arp4072_x4(silence, 800.0, 1.0)
    steady = output[len(output) // 2 :]
    assert np.all(np.isfinite(output))
    assert 0.1 < np.max(np.abs(steady)) < 13.6
    assert np.max(np.abs(np.diff(steady))) < 1.0


def test_high_drive_harmonics_match_continuous_dop853_reference():
    frequency = 1_000.0
    amplitude = 50.0
    cutoff = 2_000.0
    resonance = 0.55
    duration = 0.08
    time, reference = render_continuous_sine(
        frequency, amplitude, cutoff, resonance, duration=duration
    )
    signal = amplitude * np.sin(2.0 * np.pi * frequency * time)
    actual = dsp.arp4072_x4(signal, cutoff, resonance)

    expected_harmonics = harmonic_amplitudes(
        reference, frequency, SAMPLE_RATE, count=9, cycles=40
    )
    actual_harmonics = harmonic_amplitudes(
        actual, frequency, SAMPLE_RATE, count=9, cycles=40
    )
    for harmonic, tolerance in (
        (1, 0.025),
        (3, 0.025),
        (5, 0.025),
        (7, 0.06),
        (9, 0.06),
    ):
        index = harmonic - 1
        assert actual_harmonics[index] == pytest.approx(
            expected_harmonics[index], rel=tolerance
        )


def test_audio_rate_cutoff_and_resonance_use_polyphase_reconstruction():
    sample_count = 4_096
    time = np.arange(sample_count) / SAMPLE_RATE
    audio = 2.0 * np.sin(2.0 * np.pi * 997.0 * time)
    cutoff = 1_800.0 * np.exp2(1.5 * np.sin(2.0 * np.pi * 3_000.0 * time))
    resonance = 0.5 + 0.45 * np.sin(2.0 * np.pi * 1_700.0 * time)

    production = dsp.arp4072_controls_x4(audio, cutoff, resonance)

    audio_x4 = dsp.resampler_upsample_x4_order7(audio)
    pitch_x4 = dsp.resampler_upsample_x4_order7(np.log2(cutoff))
    resonance_x4 = np.clip(dsp.resampler_upsample_x4_order7(resonance), 0.0, 1.0)
    internal = dsp.arp4072_controls_x1(
        audio_x4,
        np.exp2(pitch_x4),
        resonance_x4,
        sample_rate=4.0 * SAMPLE_RATE,
    )
    reference = dsp.resampler_downsample_x4_order7(internal)

    np.testing.assert_allclose(production[256:], reference[256:], atol=2.0e-6)
