import numpy as np

import _triggerfish_dsp as dsp
from reference_arp4019 import circuit_values, exponential_gain

SAMPLE_RATE = 48_000.0


def _steady_rms(signal):
    steady = signal[len(signal) // 2 :]
    return np.sqrt(np.mean(steady**2))


def test_circuit_values_match_independent_component_reference():
    expected = circuit_values()
    actual = dsp.arp4019_circuit_values()
    assert actual.keys() == expected.keys()
    for name, value in expected.items():
        assert np.isclose(actual[name], value, rtol=1.0e-13, atol=1.0e-15), name


def test_linear_and_exponential_control_calibration():
    count = int(SAMPLE_RATE)
    time = np.arange(count) / SAMPLE_RATE
    audio = 1.0e-3 * np.sin(2.0 * np.pi * 200.0 * time)
    closed = np.full(count, -10.0)

    linear_unity = dsp.arp4019_x1(audio, np.full(count, 10.0), closed)
    linear_half = dsp.arp4019_x1(audio, np.full(count, 5.0), closed)
    exponential_unity = dsp.arp4019_x1(audio, closed, np.full(count, 10.0))
    exponential_minus_10db = dsp.arp4019_x1(audio, closed, np.full(count, 9.0))

    assert np.isclose(
        _steady_rms(linear_half) / _steady_rms(linear_unity), 0.5, rtol=2.0e-3
    )
    assert np.isclose(
        _steady_rms(exponential_minus_10db) / _steady_rms(exponential_unity),
        exponential_gain(9.0),
        rtol=2.0e-3,
    )
    assert np.isclose(
        _steady_rms(exponential_unity) / _steady_rms(linear_unity), 1.0, rtol=2.0e-3
    )


def test_audio_pair_compression_is_smooth_and_level_dependent():
    count = int(SAMPLE_RATE)
    time = np.arange(count) / SAMPLE_RATE
    gains = []
    maximum_steps = []
    for amplitude in (0.1, 1.0, 5.0, 10.0):
        audio = amplitude * np.sin(2.0 * np.pi * 200.0 * time)
        output = dsp.arp4019_x4(audio, np.full(count, 10.0), np.full(count, -10.0))
        steady = output[count // 2 :]
        gains.append(np.sqrt(np.mean(steady**2)) / (amplitude / np.sqrt(2.0)))
        maximum_steps.append(np.max(np.abs(np.diff(steady))))
        assert np.all(np.isfinite(output))

    assert gains[0] > gains[1] > gains[2] > gains[3]
    assert gains[1] > 0.98
    assert maximum_steps[-1] < 0.5


def test_control_current_saturates_without_a_hard_output_boundary():
    count = int(SAMPLE_RATE // 2)
    time = np.arange(count) / SAMPLE_RATE
    # A low-frequency carrier separates transfer-curve continuity from the
    # deliberately generated high harmonics at extreme control current.
    audio = 5.0 * np.sin(2.0 * np.pi * 100.0 * time)
    linear_cv = np.linspace(0.0, 100.0, count)
    exponential_cv = np.full(count, -10.0)
    output = dsp.arp4019_x4(audio, linear_cv, exponential_cv)

    assert np.all(np.isfinite(output))
    assert np.max(np.abs(output)) < 13.6
    assert np.max(np.abs(np.diff(output))) < 2.0
