import numpy as np
import pytest

import _triggerfish_dsp as dsp
from reference_diode_ladder import (
    RACK_OUTPUT_SCALE,
    STOCK_INPUT_SCALE,
    render_nonlinear_reference,
    resonance_makeup,
    transfer,
)

SAMPLE_RATE = 48_000.0


def test_stock_input_scale_matches_service_schematic_level():
    stock_saw_peak_to_peak = 6.5
    divider = 10.0 / (220.0 + 10.0)
    twice_thermal_voltage = 2.0 * 0.02585
    rack_peak_to_peak = 10.0
    expected = (
        stock_saw_peak_to_peak * divider / twice_thermal_voltage / rack_peak_to_peak
    )
    assert STOCK_INPUT_SCALE == pytest.approx(expected, rel=0.005)


def render_sine(
    processor,
    frequency,
    cutoff,
    resonance=0.0,
    high_resonance=False,
    drive_gain=0.01,
    bass=0.0,
    duration=1.5,
    sample_rate=SAMPLE_RATE,
):
    size = int(sample_rate * duration)
    time = np.arange(size) / sample_rate
    input_signal = np.sin(2.0 * np.pi * frequency * time)
    output = processor(
        input_signal,
        cutoff,
        resonance,
        high_resonance,
        drive_gain,
        bass,
        sample_rate,
    )
    return time, input_signal, output


def complex_gain(time, input_signal, output, frequency, settle=0.75):
    selected = time >= settle
    carrier = np.exp(-2j * np.pi * frequency * time[selected])
    input_phasor = np.sum(input_signal[selected] * carrier)
    output_phasor = np.sum(output[selected] * carrier)
    return output_phasor / input_phasor


def steady_rms(signal, settle_seconds=1.0):
    return np.sqrt(np.mean(signal[int(SAMPLE_RATE * settle_seconds) :] ** 2))


def saw(frequency, duration=2.0):
    time = np.arange(int(SAMPLE_RATE * duration)) / SAMPLE_RATE
    return 5.0 * (2.0 * np.mod(frequency * time, 1.0) - 1.0)


@pytest.mark.parametrize("processor", (dsp.diode_ladder_x2, dsp.diode_ladder_x4))
@pytest.mark.parametrize("frequency", (30.0, 110.0, 1_000.0, 6_000.0))
def test_small_signal_response_matches_complete_analog_model(processor, frequency):
    cutoff = 1_000.0
    time, input_signal, output = render_sine(
        processor, frequency, cutoff, resonance=0.45
    )
    measured = complex_gain(time, input_signal, output, frequency)
    # Input drive and output calibration are deliberately independent. Their
    # product scales the raw circuit response in the linear limit; resonance
    # makeup is then applied as a separate musical output-stage decision.
    expected = (
        STOCK_INPUT_SCALE
        * RACK_OUTPUT_SCALE
        * 0.01
        * resonance_makeup(0.45)
        * transfer(frequency, cutoff, resonance=0.45)
    )
    assert abs(measured) == pytest.approx(abs(expected), rel=0.035, abs=2.0e-6)


@pytest.mark.parametrize("sample_rate", (44_100.0, 96_000.0, 192_000.0))
def test_small_signal_response_is_consistent_across_host_rates(sample_rate):
    frequency = 1_000.0
    cutoff = 1_000.0
    time, input_signal, output = render_sine(
        dsp.diode_ladder_x4,
        frequency,
        cutoff,
        resonance=0.45,
        duration=1.0,
        sample_rate=sample_rate,
    )
    measured = complex_gain(time, input_signal, output, frequency, settle=0.5)
    expected = (
        STOCK_INPUT_SCALE
        * RACK_OUTPUT_SCALE
        * 0.01
        * resonance_makeup(0.45)
        * transfer(frequency, cutoff, resonance=0.45)
    )
    assert abs(measured) == pytest.approx(abs(expected), rel=0.035, abs=2.0e-6)


@pytest.mark.parametrize("processor", (dsp.diode_ladder_x2, dsp.diode_ladder_x4))
def test_bass_control_reaches_documented_four_db_extension(processor):
    frequency = 32.0
    time, input_signal, stock = render_sine(processor, frequency, 1_000.0, bass=0.0)
    _, _, full = render_sine(processor, frequency, 1_000.0, bass=1.0)
    stock_gain = abs(complex_gain(time, input_signal, stock, frequency))
    full_gain = abs(complex_gain(time, input_signal, full, frequency))
    extension_db = 20.0 * np.log10(full_gain / stock_gain)
    assert extension_db == pytest.approx(4.0, abs=0.2)


@pytest.mark.parametrize("frequency", (55.0, 110.0, 220.0, 440.0))
def test_stock_open_filter_keeps_nominal_rack_oscillator_level(frequency):
    signal = saw(frequency)
    output = dsp.diode_ladder_x4(signal, 20_000.0, 0.0, False, 1.0, 0.0, SAMPLE_RATE)
    selected = output[int(SAMPLE_RATE) :]
    level_ratio = steady_rms(output) / steady_rms(signal)

    # Circuit coupling and nonlinear compression make exact unity neither
    # possible nor desirable, but a normal +/-5 V oscillator must stay in the
    # same practical Rack voltage range when the filter is open.
    assert 0.7 < level_ratio < 1.2
    assert 4.0 < np.max(np.abs(selected)) < 7.0


def test_resonance_makeup_retains_driven_saw_without_erasing_authentic_thinning():
    signal = saw(110.0)
    plain = dsp.diode_ladder_x4(signal, 1_000.0, 0.0, False, 1.0, 0.0, SAMPLE_RATE)
    stock = dsp.diode_ladder_x4(signal, 1_000.0, 1.0, False, 1.0, 0.0, SAMPLE_RATE)
    high = dsp.diode_ladder_x4(signal, 1_000.0, 1.0, True, 1.0, 0.0, SAMPLE_RATE)
    stock_ratio = steady_rms(stock) / steady_rms(plain)
    high_ratio = steady_rms(high) / steady_rms(plain)

    assert 0.5 < stock_ratio < 0.8
    assert 0.4 < high_ratio < stock_ratio


@pytest.mark.parametrize("processor", (dsp.diode_ladder_x2, dsp.diode_ladder_x4))
def test_extreme_modulation_stays_finite(processor):
    size = 32_768
    time = np.arange(size) / SAMPLE_RATE
    signal = 10.0 * (
        np.sin(2.0 * np.pi * 997.0 * time) + 0.5 * np.sin(2.0 * np.pi * 7_013.0 * time)
    )
    output = processor(signal, 19_000.0, 1.0, True, 66.6, 1.0, SAMPLE_RATE)
    assert np.isfinite(output).all()
    assert np.max(np.abs(output)) <= 20.0001


def test_stock_resonance_decays_and_high_mode_reaches_documented_threshold():
    size = int(3.0 * SAMPLE_RATE)
    impulse = np.zeros(size)
    impulse[0] = 5.0
    stock = dsp.diode_ladder_x4(impulse, 1_000.0, 1.0, False, 1.0, 0.0, SAMPLE_RATE)
    high_below = dsp.diode_ladder_x4(
        impulse, 1_000.0, 0.67, True, 1.0, 0.0, SAMPLE_RATE
    )
    high_above = dsp.diode_ladder_x4(
        impulse, 1_000.0, 0.70, True, 1.0, 0.0, SAMPLE_RATE
    )
    stock_tail = np.sqrt(np.mean(stock[-int(SAMPLE_RATE) :] ** 2))
    below_tail = np.sqrt(np.mean(high_below[-int(SAMPLE_RATE) :] ** 2))
    above_tail = np.sqrt(np.mean(high_above[-int(SAMPLE_RATE) :] ** 2))
    assert stock_tail < 1.0e-5
    assert below_tail < 1.0e-5
    assert above_tail > 0.05


def test_complete_coupling_network_rejects_dc():
    signal = np.full(int(4.0 * SAMPLE_RATE), 5.0)
    output = dsp.diode_ladder_x2(signal, 1_000.0, 0.5, False, 1.0, 0.0, SAMPLE_RATE)
    assert np.max(np.abs(output[-1_000:])) < 1.0e-5


def test_nonlinear_harmonics_match_independent_continuous_reference():
    frequency = 220.0
    duration = 0.2
    size = int(duration * SAMPLE_RATE)
    time = np.arange(size) / SAMPLE_RATE
    input_signal = 5.0 * np.sin(2.0 * np.pi * frequency * time)
    production = dsp.diode_ladder_x4(
        input_signal, 1_500.0, 0.55, False, 10.0, 0.5, SAMPLE_RATE
    )
    reference_time, reference = render_nonlinear_reference(
        lambda value: 5.0 * np.sin(2.0 * np.pi * frequency * value),
        SAMPLE_RATE,
        duration,
        1_500.0,
        resonance=0.55,
        drive_gain=10.0,
        bass=0.5,
    )
    assert np.array_equal(reference_time, time)

    # The 0.1 s analysis window contains exactly 22 periods, preventing the
    # dominant fundamental from leaking into the much smaller harmonics.
    selected = time >= 0.1
    for harmonic in (1, 3, 5, 7, 9):
        carrier = np.exp(-2j * np.pi * harmonic * frequency * time[selected])
        production_level = abs(np.sum(production[selected] * carrier))
        reference_level = abs(np.sum(reference[selected] * carrier))
        assert production_level == pytest.approx(reference_level, rel=0.035, abs=2.0e-3)


def test_two_and_four_times_outputs_agree_at_normal_settings():
    size = 48_000
    time = np.arange(size) / SAMPLE_RATE
    signal = 3.0 * np.sin(2.0 * np.pi * 220.0 * time)
    x2 = dsp.diode_ladder_x2(signal, 1_200.0, 0.6, False, 1.0, 0.5, SAMPLE_RATE)
    x4 = dsp.diode_ladder_x4(signal, 1_200.0, 0.6, False, 1.0, 0.5, SAMPLE_RATE)
    x2_gain = complex_gain(time, signal, x2, 220.0, settle=0.25)
    x4_gain = complex_gain(time, signal, x4, 220.0, settle=0.25)
    # The cascaded 4x resampler has approximately one more host sample of
    # group delay. Compare steady-state magnitude rather than unaligned samples.
    assert abs(x2_gain) == pytest.approx(abs(x4_gain), rel=0.002)


def test_vca_runs_inside_oversampled_path_without_changing_nominal_level():
    size = int(1.5 * SAMPLE_RATE)
    time = np.arange(size) / SAMPLE_RATE
    signal = 5.0 * np.sin(2.0 * np.pi * 1_000.0 * time)
    control = np.ones(size)
    host_rate = dsp.diode_ladder_vca_x4(
        signal, control, 20_000.0, 0.0, False, 1.0, 0.0, SAMPLE_RATE, False
    )
    oversampled = dsp.diode_ladder_vca_x4(
        signal, control, 20_000.0, 0.0, False, 1.0, 0.0, SAMPLE_RATE, True
    )
    selected = slice(int(SAMPLE_RATE), None)
    host_level = np.sqrt(np.mean(host_rate[selected] ** 2))
    oversampled_level = np.sqrt(np.mean(oversampled[selected] ** 2))
    assert np.isfinite(oversampled).all()
    assert np.max(np.abs(oversampled)) <= 12.0001
    assert oversampled_level == pytest.approx(host_level, rel=0.002)
