import numpy as np
import pytest
from scipy.signal import resample_poly

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
    stock_saw_peak_to_peak = 5.5
    divider = 2.2 / (220.0 + 2.2)
    twice_thermal_voltage = 2.0 * 0.02585
    rack_peak_to_peak = 10.0
    expected = (
        stock_saw_peak_to_peak * divider / twice_thermal_voltage / rack_peak_to_peak
    )
    assert STOCK_INPUT_SCALE == pytest.approx(expected, rel=0.005)


def test_stock_saw_stays_near_the_ladder_small_signal_region():
    rack_saw_peak = 5.0
    normalized_peak = rack_saw_peak * STOCK_INPUT_SCALE

    # A stock saw reaches tanh(0.527), leaving deliberate junction curvature
    # without starting close to the asymptotic rail.
    assert normalized_peak == pytest.approx(0.527, rel=0.005)
    assert np.tanh(normalized_peak) == pytest.approx(0.483, rel=0.005)


def test_linear_fm_cutoff_pinches_off_smoothly():
    requested = np.asarray((-20.0, -5.0, 0.0, 5.0, 20.0, 19_990.0, 20_000.0))
    mapped = np.asarray([dsp.diode_ladder_map_cutoff(value) for value in requested])

    assert np.all(np.diff(mapped) > 0.0)
    assert mapped[0] > 0.0
    assert mapped[2] == pytest.approx(np.log(2.0), rel=1.0e-5)
    assert mapped[4] == pytest.approx(20.0, abs=1.0e-6)
    assert mapped[-1] < 20_000.0


def test_stock_oscillator_does_not_start_in_deep_ladder_saturation():
    duration = 2.0
    size = int(SAMPLE_RATE * duration)
    pitch = np.full(size, np.log2(82.4069 / 261.625565))
    zero = np.zeros(size)
    oscillator = dsp.tb303_oscillator_x2(pitch, zero, zero, zero, zero, SAMPLE_RATE)[
        :, 2
    ]
    small = dsp.diode_ladder_x4(
        0.01 * oscillator, 92.5, 0.0, False, 1.0, 0.0, SAMPLE_RATE
    )
    stock = dsp.diode_ladder_x4(oscillator, 92.5, 0.0, False, 1.0, 0.0, SAMPLE_RATE)
    selected = slice(int(SAMPLE_RATE), None)
    linear = small[selected] / 0.01
    measured = stock[selected]
    fitted_gain = np.dot(measured, linear) / np.dot(linear, linear)
    residual = measured - fitted_gain * linear
    residual_db = 20.0 * np.log10(
        np.sqrt(np.mean(residual**2)) / np.sqrt(np.mean(measured**2))
    )

    assert residual_db < -24.0


def test_low_cutoff_suppresses_saw_reset_resampler_ringing():
    duration = 1.0
    size = int(SAMPLE_RATE * duration)
    frequency = 82.4069
    pitch = np.full(size, np.log2(frequency / 261.625565))
    zero = np.zeros(size)
    candidate = dsp.tb303_oscillator_x2(pitch, zero, zero, zero, zero, SAMPLE_RATE)[
        :, 2
    ]

    reference_size = 16 * size
    reference_pitch = np.full(reference_size, np.log2(frequency / 261.625565))
    reference_zero = np.zeros(reference_size)
    high_rate = dsp.tb303_oscillator_x1(
        reference_pitch,
        reference_zero,
        reference_zero,
        reference_zero,
        reference_zero,
        16.0 * SAMPLE_RATE,
    )[:, 2]
    reference = resample_poly(high_rate, 1, 16, window=("kaiser", 12.0))[:size]

    candidate_filtered = dsp.diode_ladder_x4(
        candidate, 92.5, 0.0, False, 1.0, 0.0, SAMPLE_RATE
    )
    reference_filtered = dsp.diode_ladder_x4(
        reference, 92.5, 0.0, False, 1.0, 0.0, SAMPLE_RATE
    )
    selected = slice(size // 2, None)
    measured = candidate_filtered[selected]
    expected = reference_filtered[selected]
    fitted_gain = np.dot(measured, expected) / np.dot(measured, measured)
    residual = fitted_gain * measured - expected
    residual_db = 20.0 * np.log10(
        np.sqrt(np.mean(residual**2)) / np.sqrt(np.mean(expected**2))
    )

    assert residual_db < -35.0


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
    assert 0.65 < level_ratio < 1.1
    assert 4.0 < np.max(np.abs(selected)) < 8.25


def test_resonance_makeup_retains_driven_saw_without_erasing_authentic_thinning():
    signal = saw(110.0)
    plain = dsp.diode_ladder_x4(signal, 1_000.0, 0.0, False, 1.0, 0.0, SAMPLE_RATE)
    stock = dsp.diode_ladder_x4(signal, 1_000.0, 1.0, False, 1.0, 0.0, SAMPLE_RATE)
    high = dsp.diode_ladder_x4(signal, 1_000.0, 1.0, True, 1.0, 0.0, SAMPLE_RATE)
    stock_ratio = steady_rms(stock) / steady_rms(plain)
    high_ratio = steady_rms(high) / steady_rms(plain)

    assert 0.35 < stock_ratio < 0.65
    assert stock_ratio < high_ratio < 1.5


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


@pytest.mark.parametrize("processor", (dsp.diode_ladder_x2, dsp.diode_ladder_x4))
def test_bass_extension_rejects_nonlinear_saw_dc_without_rail_clipping(processor):
    signal = saw(110.0)
    output = processor(signal, 500.0, 0.4, False, 66.6, 1.0, SAMPLE_RATE)
    settled = output[int(SAMPLE_RATE) :]

    assert abs(np.mean(settled)) < 1.0e-3
    assert np.max(np.abs(settled)) < 10.0
    assert not np.any(np.abs(settled) == 20.0)


def test_drive_and_bass_moves_do_not_create_discontinuous_clicks():
    duration = 3.0
    time = np.arange(int(SAMPLE_RATE * duration)) / SAMPLE_RATE
    signal = 5.0 * (2.0 * np.mod(110.0 * time, 1.0) - 1.0)
    drive = np.where(time < 1.0, 1.0, 16.0)
    bass = np.where(time < 1.0, 0.0, 1.0)
    output = dsp.diode_ladder_modulated_x4(
        signal, drive, bass, 500.0, 0.4, False, SAMPLE_RATE
    )
    delta = np.abs(np.diff(output))
    transition = int(SAMPLE_RATE)

    baseline = np.max(delta[transition - 5_000 : transition - 100])
    assert delta[transition - 1] < baseline
    assert np.max(delta[transition - 5 : transition + 10]) < 0.01
    assert abs(np.mean(output[int(2.0 * SAMPLE_RATE) :])) < 1.0e-3


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
        assert production_level == pytest.approx(reference_level, rel=0.065, abs=2.0e-3)


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


def test_low_octave_bass_and_resonance_overload_has_no_hard_rail():
    size = int(4.0 * SAMPLE_RATE)
    pitch = np.full(size, np.log2(10.3 / 261.625565))
    zero = np.zeros(size)
    oscillator = dsp.tb303_oscillator_x2(pitch, zero, zero, zero, zero, SAMPLE_RATE)[
        :, 2
    ]
    control = np.ones(size)
    output = dsp.diode_ladder_vca_x4(
        oscillator,
        control,
        2_000.0,
        0.589,
        True,
        1.0,
        0.472,
        SAMPLE_RATE,
        True,
    )
    settled = output[int(2.0 * SAMPLE_RATE) :]

    assert np.isfinite(settled).all()
    assert 8.5 < np.max(np.abs(settled)) < 11.1
    assert not np.any(np.abs(settled) == 12.0)
    assert np.count_nonzero(np.diff(settled) == 0.0) < 0.001 * settled.size


def test_extreme_modulated_voice_uses_smooth_safety_compliance():
    size = 32_768
    rng = np.random.default_rng(303)
    audio = rng.choice((-10.0, 10.0), size)
    cutoff = 10.0 ** rng.uniform(np.log10(5.0), np.log10(20_000.0), size)
    base_control = rng.integers(0, 2, size).astype(float)
    accent_control = rng.integers(0, 2, size).astype(float)
    rendered = dsp.diode_ladder_voice_x4(
        audio,
        cutoff,
        base_control,
        accent_control,
        1.0,
        True,
        66.6,
        1.0,
        SAMPLE_RATE,
    )

    assert rendered[-1, 2] == 0.0
    assert np.isfinite(rendered[:, :2]).all()
    assert 11.0 < np.max(np.abs(rendered[:, 1])) < 12.0
    assert not np.any(np.abs(rendered[:, :2]) == 12.0)


def test_two_and_four_times_vca_overload_remain_consistent():
    size = int(2.0 * SAMPLE_RATE)
    time = np.arange(size) / SAMPLE_RATE
    signal = 10.0 * np.sin(2.0 * np.pi * 997.0 * time)
    control = np.ones(size)
    arguments = (signal, control, 15_000.0, 0.8, True, 66.6, 0.0, SAMPLE_RATE)
    x2 = dsp.diode_ladder_vca_x2(*arguments, True)
    x4 = dsp.diode_ladder_vca_x4(*arguments, True)
    selected = slice(int(SAMPLE_RATE), None)
    x2_rms = np.sqrt(np.mean(x2[selected] ** 2))
    x4_rms = np.sqrt(np.mean(x4[selected] ** 2))

    assert np.isfinite(x2).all()
    assert np.isfinite(x4).all()
    assert x2_rms == pytest.approx(x4_rms, rel=0.015)

    def off_harmonic_level(output):
        # The selected one-second window contains an integer number of 997 Hz
        # cycles. Remove the physically valid in-band odd harmonics; remaining
        # FFT energy is dominated by folded harmonics and numerical noise.
        spectrum = 2.0 * np.abs(np.fft.rfft(output[selected])) / int(SAMPLE_RATE)
        expected = np.zeros(spectrum.size, dtype=bool)
        expected[:2] = True
        for harmonic in range(1, 25, 2):
            bin_index = int(harmonic * 997)
            if bin_index < spectrum.size:
                expected[bin_index - 1 : bin_index + 2] = True
        residual = np.sqrt(np.sum(spectrum[~expected] ** 2))
        return 20.0 * np.log10(residual / spectrum[997])

    x2_residual = off_harmonic_level(x2)
    x4_residual = off_harmonic_level(x4)
    assert x4_residual < -45.0
    assert x4_residual < x2_residual - 8.0


def test_audio_rate_pitch_fm_and_resonance_use_polyphase_reconstruction():
    sample_count = 4_096
    time = np.arange(sample_count) / SAMPLE_RATE
    audio = 2.0 * np.sin(2.0 * np.pi * 997.0 * time)
    cutoff = 1_200.0 * np.exp2(0.7 * np.sin(2.0 * np.pi * 2_200.0 * time))
    linear_fm = 900.0 * np.sin(2.0 * np.pi * 3_100.0 * time)
    resonance = 0.45 + 0.4 * np.sin(2.0 * np.pi * 1_300.0 * time)

    production = dsp.diode_ladder_controls_x4(audio, cutoff, linear_fm, resonance)

    audio_x4 = dsp.resampler_upsample_x4_order7(audio)
    pitch_x4 = dsp.resampler_upsample_x4_order7(np.log2(cutoff))
    linear_fm_x4 = dsp.resampler_upsample_x4_order7(linear_fm)
    resonance_x4 = np.clip(dsp.resampler_upsample_x4_order7(resonance), 0.0, 1.0)
    internal = dsp.diode_ladder_controls_x1(
        audio_x4,
        np.exp2(pitch_x4),
        linear_fm_x4,
        resonance_x4,
        sample_rate=4.0 * SAMPLE_RATE,
    )
    reference = dsp.resampler_downsample_x4_order7(internal)

    np.testing.assert_allclose(production[256:], reference[256:], atol=2.0e-6)
