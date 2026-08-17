import numpy as np
import pytest
from scipy.integrate import quad
from scipy.signal import resample_poly
from scipy.signal.windows import hann
from scipy.special import lambertw

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000
REFERENCE_FACTOR = 16
CHARACTERS = (0, 1, 2)


def lockhart_reference(inputs):
    thermal_voltage = 0.025864
    resistance = 15_000.0
    load_resistance = 7_500.0
    delta = load_resistance * 1.0e-17 / thermal_voltage
    beta = (2.0 * load_resistance + resistance) / (thermal_voltage * resistance)
    value = np.asarray(inputs) * 1.069465474121151 / 3.0
    for _ in range(4):
        polarity = np.where(value >= 0.0, 1.0, -1.0)
        value -= (
            polarity
            * thermal_voltage
            * lambertw(delta * np.exp(polarity * beta * value)).real
        )
    return 3.0 * 0.9350465482253767 * value


def serge_reference(inputs):
    junction_voltage = 1.752 * 0.025864
    series_drop = 33_000.0 * 2.52e-9
    delta = series_drop / junction_voltage
    value = np.asarray(inputs) * 0.33073419047751157
    for _ in range(6):
        polarity = np.where(value >= 0.0, 1.0, -1.0)
        w = lambertw(
            delta * np.exp((polarity * value + series_drop) / junction_voltage)
        ).real
        value -= 2.0 * polarity * (junction_voltage * w - series_drop)
    return 4.0 * 0.7727253519539319 * value


def constant_render(
    function, frequency, morph, fold, symmetry, seconds, *, adaa=True, character=0
):
    samples = round(SAMPLE_RATE * seconds)
    frequency_signal = np.full(samples, frequency, dtype=np.float64)
    return function(
        frequency_signal,
        np.full(samples, morph),
        np.full(samples, fold),
        np.full(samples, symmetry),
        sample_rate=SAMPLE_RATE,
        adaa=adaa,
        character=character,
    )


def spectral_error(candidate, reference):
    window = np.hanning(reference.size)
    candidate_spectrum = np.abs(np.fft.rfft(candidate * window))
    reference_spectrum = np.abs(np.fft.rfft(reference * window))
    candidate_spectrum /= np.linalg.norm(candidate_spectrum)
    reference_spectrum /= np.linalg.norm(reference_spectrum)
    return np.linalg.norm(candidate_spectrum - reference_spectrum)


def render_unison(voices, spread, seconds=2.0):
    positions = dsp.unison_pitch_positions(voices)
    spread_cents = dsp.unison_spread_cents(spread)
    rendered = [
        constant_render(
            dsp.wavefold_oscillator_x4,
            440.0 * 2.0 ** (spread_cents * position / 1200.0),
            0.5,
            0.6,
            0.0,
            seconds,
            adaa=False,
            character=1,
        )
        for position in positions
    ]
    return dsp.unison_output_gain(voices) * np.sum(rendered, axis=0)


def test_wavefolder_tables_are_odd_and_their_primitives_are_consistent():
    inputs = np.linspace(-12.0, 12.0, 24_001)
    for character in CHARACTERS:
        transfer = dsp.wavefolder_transfer(inputs, character)
        primitive = dsp.wavefolder_primitive(inputs, character)

        np.testing.assert_allclose(transfer, -transfer[::-1], atol=2e-12)
        numerical_derivative = np.gradient(primitive, inputs)
        np.testing.assert_allclose(
            numerical_derivative[2:-2], transfer[2:-2], atol=2e-5
        )


@pytest.mark.parametrize(
    ("voices", "expected"),
    [
        (1, [0.0]),
        (2, [-1.0, 1.0]),
        (3, [-1.0, 0.0, 1.0]),
        (4, [-1.0, -0.18, 0.20, 0.98]),
    ],
)
def test_unison_pitch_layouts_are_centred_and_nonuniform(voices, expected):
    positions = dsp.unison_pitch_positions(voices)
    np.testing.assert_allclose(positions, expected, atol=1e-15)
    assert np.mean(positions) == pytest.approx(0.0, abs=1e-15)


def test_unison_energy_normalization_keeps_detuned_output_in_the_same_range():
    reference = render_unison(1, 0.5)
    reference_rms = np.sqrt(np.mean(reference[SAMPLE_RATE // 2 :] ** 2))
    for voices in (2, 3, 4):
        output = render_unison(voices, 0.5)
        rms = np.sqrt(np.mean(output[SAMPLE_RATE // 2 :] ** 2))
        assert 0.8 < rms / reference_rms < 1.2


def test_circuit_tables_follow_the_lockhart_and_serge_equations():
    inputs = np.linspace(-10.0, 10.0, 20_001)
    np.testing.assert_allclose(
        dsp.wavefolder_transfer(inputs, 1), lockhart_reference(inputs), atol=3e-8
    )
    np.testing.assert_allclose(
        dsp.wavefolder_transfer(inputs, 2), serge_reference(inputs), atol=2e-7
    )


def test_first_order_adaa_is_the_segment_average_of_the_transfer_curve():
    inputs = np.array(
        [0.1, 1.7, -2.3, 4.9, -7.1, -7.1, 0.000001, 3.2, -0.4],
        dtype=np.float64,
    )
    for character in CHARACTERS:
        output = dsp.wavefolder_adaa(inputs, character)
        assert output[0] == dsp.wavefolder_transfer(inputs[:1], character)[0]
        for index in range(1, inputs.size):
            left = inputs[index - 1]
            right = inputs[index]
            if left == right:
                expected = dsp.wavefolder_transfer(np.array([left]), character)[0]
            else:
                integral, _ = quad(
                    lambda value: dsp.wavefolder_transfer(np.array([value]), character)[
                        0
                    ],
                    left,
                    right,
                    epsabs=1e-9,
                    epsrel=1e-9,
                    limit=500,
                )
                expected = integral / (right - left)
            assert output[index] == pytest.approx(expected, abs=2e-9)


def test_wavefold_oscillator_is_finite_with_audio_rate_control_modulation():
    samples = SAMPLE_RATE
    time = np.arange(samples) / SAMPLE_RATE
    frequency = 880.0 * np.exp2(0.9 * np.sin(2.0 * np.pi * 997.0 * time))
    morph = 0.5 + 0.5 * np.sin(2.0 * np.pi * 733.0 * time + 0.2)
    fold = 0.5 + 0.5 * np.sin(2.0 * np.pi * 601.0 * time + 0.4)
    symmetry = 0.8 * np.sin(2.0 * np.pi * 431.0 * time + 0.6)

    for character in CHARACTERS:
        for function in (dsp.wavefold_oscillator_x2, dsp.wavefold_oscillator_x4):
            output = function(
                frequency, morph, fold, symmetry, SAMPLE_RATE, True, character
            )
            assert np.all(np.isfinite(output))
            assert np.max(np.abs(output)) < 4.0


def test_zero_fold_is_an_exact_character_independent_bypass():
    samples = SAMPLE_RATE
    frequency = np.full(samples, 251.0)
    fold = np.zeros(samples)
    symmetry = np.zeros(samples)
    for morph in (0.0, 1.0):
        outputs = [
            dsp.wavefold_oscillator_x4(
                frequency,
                np.full(samples, morph),
                fold,
                symmetry,
                SAMPLE_RATE,
                False,
                character,
            )
            for character in CHARACTERS
        ]
        np.testing.assert_array_equal(outputs[1], outputs[0])
        np.testing.assert_array_equal(outputs[2], outputs[0])


def test_unfolded_triangle_has_the_ideal_low_harmonic_spectrum():
    frequency = 250.0
    output = constant_render(
        dsp.wavefold_oscillator_x4,
        frequency,
        1.0,
        0.0,
        0.0,
        2.0,
        adaa=False,
        character=1,
    )[SAMPLE_RATE:]
    spectrum = np.abs(np.fft.rfft(output))
    fundamental = spectrum[round(frequency)]
    for harmonic in (3, 5, 7, 9):
        measured = spectrum[round(harmonic * frequency)] / fundamental
        assert measured == pytest.approx(1.0 / harmonic**2, rel=0.015)


def test_x4_improves_the_high_note_fold_spectrum_over_x2():
    frequency = 4_186.0
    padding = 0.125
    duration = 0.5
    total = padding + duration
    high_rate = SAMPLE_RATE * REFERENCE_FACTOR
    high_samples = round(total * high_rate)
    start = round(padding * SAMPLE_RATE)
    stop = start + round(duration * SAMPLE_RATE)
    for character in CHARACTERS:
        reference = dsp.wavefold_oscillator_x1(
            np.full(high_samples, frequency),
            np.full(high_samples, 0.5),
            np.ones(high_samples),
            np.zeros(high_samples),
            high_rate,
            False,
            character,
        )
        reference = resample_poly(
            reference, 1, REFERENCE_FACTOR, window=("kaiser", 12.0)
        )[start:stop]
        x2 = constant_render(
            dsp.wavefold_oscillator_x2,
            frequency,
            0.5,
            1.0,
            0.0,
            total,
            character=character,
        )[start:stop]
        x4 = constant_render(
            dsp.wavefold_oscillator_x4,
            frequency,
            0.5,
            1.0,
            0.0,
            total,
            character=character,
        )[start:stop]

        assert spectral_error(x4, reference) < spectral_error(x2, reference)


def off_harmonic_energy_db(signal, frequency):
    spectrum = np.abs(np.fft.rfft(signal * hann(signal.size, sym=False))) ** 2
    harmonic_bins = np.zeros(spectrum.size, dtype=bool)
    harmonic_bins[:3] = True
    for harmonic in range(1, int((SAMPLE_RATE / 2) // frequency) + 1):
        center = round(harmonic * frequency)
        harmonic_bins[max(0, center - 2) : center + 3] = True
    ratio = np.sum(spectrum[~harmonic_bins]) / np.sum(spectrum)
    return 10.0 * np.log10(max(ratio, np.finfo(float).tiny))


def test_adaa_and_x16_reduce_high_note_external_folder_aliases():
    frequency = 4_187.0
    settling_samples = SAMPLE_RATE // 2
    analysis_samples = SAMPLE_RATE
    samples = settling_samples + analysis_samples
    time = np.arange(samples) / SAMPLE_RATE
    audio = np.sin(2.0 * np.pi * frequency * time)
    fold = np.ones(samples)
    symmetry = np.zeros(samples)

    for character in CHARACTERS:
        direct_x4 = dsp.wavefolder_external_x4(
            audio, fold, symmetry, SAMPLE_RATE, False, character
        )[settling_samples:]
        adaa_x4 = dsp.wavefolder_external_x4(
            audio, fold, symmetry, SAMPLE_RATE, True, character
        )[settling_samples:]
        direct_x16 = dsp.wavefolder_external_x16(
            audio, fold, symmetry, SAMPLE_RATE, False, character
        )[settling_samples:]
        direct_alias = off_harmonic_energy_db(direct_x4, frequency)
        assert off_harmonic_energy_db(adaa_x4, frequency) < direct_alias - 2.0
        assert off_harmonic_energy_db(direct_x16, frequency) < direct_alias - 2.0
