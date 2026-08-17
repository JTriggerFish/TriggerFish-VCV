"""Compare wavefolder oversampling quality and CPU cost."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import time

import numpy as np
from scipy.signal import resample_poly
from scipy.signal.windows import hann

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000
REFERENCE_FACTOR = 16
REFERENCE_WINDOW = ("kaiser", 12.0)
VARIANTS = {
    "2x direct": (dsp.wavefold_oscillator_x2, False),
    "4x direct": (dsp.wavefold_oscillator_x4, False),
    "4x + ADAA": (dsp.wavefold_oscillator_x4, True),
    "16x direct": (dsp.wavefold_oscillator_x16, False),
}
CHARACTERS = {"Hinge": 0, "Lockhart": 1, "Serge": 2}


@dataclass(frozen=True)
class Scenario:
    name: str
    frequency: float
    morph: float
    fold: float
    symmetry: float = 0.0
    modulation: str | None = None


SCENARIOS = (
    Scenario("C4, medium fold", 261.63, 0.5, 0.5),
    Scenario("C6, full fold", 1_046.5, 1.0, 1.0),
    Scenario("C7, full fold", 2_093.0, 0.5, 1.0),
    Scenario("C8, full fold", 4_186.0, 0.5, 1.0),
    Scenario("C6, asymmetric full fold", 1_046.5, 0.25, 1.0, 0.6),
    Scenario("997 Hz exponential FM", 880.0, 0.5, 0.8, 0.3, "pitch"),
    Scenario("733 Hz waveform morph", 880.0, 0.5, 0.8, 0.3, "morph"),
    Scenario("601 Hz fold modulation", 880.0, 0.5, 0.8, 0.3, "fold"),
    Scenario("431 Hz symmetry modulation", 880.0, 0.5, 0.8, 0.3, "symmetry"),
)


def controls(scenario: Scenario, sample_rate: int, seconds: float):
    samples = round(sample_rate * seconds)
    time_axis = np.arange(samples, dtype=np.float64) / sample_rate
    if scenario.modulation is None:
        return tuple(
            np.full(samples, value, dtype=np.float64)
            for value in (
                scenario.frequency,
                scenario.morph,
                scenario.fold,
                scenario.symmetry,
            )
        )
    frequency = np.full(samples, scenario.frequency, dtype=np.float64)
    morph = np.full(samples, scenario.morph, dtype=np.float64)
    fold = np.full(samples, scenario.fold, dtype=np.float64)
    symmetry = np.full(samples, scenario.symmetry, dtype=np.float64)
    if scenario.modulation == "pitch":
        frequency *= np.exp2(0.7 * np.sin(2.0 * np.pi * 997.0 * time_axis))
    elif scenario.modulation == "morph":
        morph += 0.45 * np.sin(2.0 * np.pi * 733.0 * time_axis + 0.2)
    elif scenario.modulation == "fold":
        fold += 0.2 * np.sin(2.0 * np.pi * 601.0 * time_axis + 0.4)
    elif scenario.modulation == "symmetry":
        symmetry += 0.5 * np.sin(2.0 * np.pi * 431.0 * time_axis + 0.6)
    return frequency, morph, fold, symmetry


def render(function, adaa, character, scenario, sample_rate, seconds):
    return function(
        *controls(scenario, sample_rate, seconds),
        sample_rate=sample_rate,
        adaa=adaa,
        character=character,
    )


def spectral_error(candidate, reference):
    """Normalized magnitude-spectrum error in dB; more negative is better."""
    window = np.hanning(reference.size)
    candidate_spectrum = np.abs(np.fft.rfft(candidate * window))
    reference_spectrum = np.abs(np.fft.rfft(reference * window))
    candidate_spectrum /= np.linalg.norm(candidate_spectrum)
    reference_spectrum /= np.linalg.norm(reference_spectrum)
    error = np.linalg.norm(candidate_spectrum - reference_spectrum)
    return 20.0 * np.log10(max(error, np.finfo(float).tiny))


def quality_results(seconds):
    padding = 0.25
    total = seconds + 2.0 * padding
    start = round(padding * SAMPLE_RATE)
    stop = start + round(seconds * SAMPLE_RATE)
    results = []
    for character_name, character in CHARACTERS.items():
        for scenario in SCENARIOS:
            reference = render(
                dsp.wavefold_oscillator_x1,
                False,
                character,
                scenario,
                SAMPLE_RATE * REFERENCE_FACTOR,
                total,
            )
            reference = resample_poly(
                reference, 1, REFERENCE_FACTOR, window=REFERENCE_WINDOW
            )[start:stop]
            for name, (function, adaa) in VARIANTS.items():
                candidate = render(
                    function, adaa, character, scenario, SAMPLE_RATE, total
                )[start:stop]
                results.append(
                    (
                        character_name,
                        scenario.name,
                        name,
                        spectral_error(candidate, reference),
                    )
                )
    return results


def cpu_result(function, adaa, character, scenario, samples, repeats):
    arguments = controls(scenario, SAMPLE_RATE, samples / SAMPLE_RATE)
    function(*arguments, SAMPLE_RATE, adaa, character)
    timings = []
    checksum = 0.0
    for _ in range(repeats):
        started = time.perf_counter_ns()
        output = function(*arguments, SAMPLE_RATE, adaa, character)
        timings.append(time.perf_counter_ns() - started)
        checksum += float(output[-1])
    if not np.isfinite(checksum):
        raise RuntimeError("non-finite benchmark output")
    return float(np.median(timings)) / samples


def timbre_results():
    frequency = 100.0
    samples = 2 * SAMPLE_RATE
    analysis = slice(SAMPLE_RATE, None)
    frequency_signal = np.full(samples, frequency)
    zero = np.zeros(samples)
    results = []
    for fold in (0.25, 0.5, 1.0):
        for character_name, character in CHARACTERS.items():
            output = dsp.wavefold_oscillator_x4(
                frequency_signal,
                zero,
                np.full(samples, fold),
                zero,
                SAMPLE_RATE,
                False,
                character,
            )[analysis]
            spectrum = np.abs(np.fft.rfft(output))
            harmonics = np.array(
                [spectrum[round(frequency * index)] for index in range(1, 31)]
            )
            harmonic_numbers = np.arange(1, harmonics.size + 1)
            centroid = (
                frequency
                * np.sum(harmonic_numbers * harmonics**2)
                / np.sum(harmonics**2)
            )
            relative = 20.0 * np.log10(
                np.maximum(harmonics / harmonics[0], np.finfo(float).tiny)
            )
            results.append(
                (
                    character_name,
                    fold,
                    np.sqrt(np.mean(output**2)),
                    centroid,
                    relative[2],
                    relative[4],
                    relative[6],
                )
            )
    return results


def off_harmonic_energy_db(signal, frequency):
    """Energy outside the continuous-time harmonics below Nyquist."""
    spectrum = np.abs(np.fft.rfft(signal * hann(signal.size, sym=False))) ** 2
    harmonic_bins = np.zeros(spectrum.size, dtype=bool)
    harmonic_bins[:3] = True
    for harmonic in range(1, int((SAMPLE_RATE / 2) // frequency) + 1):
        center = round(harmonic * frequency)
        harmonic_bins[max(0, center - 2) : center + 3] = True
    ratio = np.sum(spectrum[~harmonic_bins]) / np.sum(spectrum)
    return 10.0 * np.log10(max(ratio, np.finfo(float).tiny))


def alias_results():
    """Stress the external folder with a coherent high sine and full fold."""
    frequency = 4_187.0
    settling_samples = SAMPLE_RATE // 2
    analysis_samples = SAMPLE_RATE
    samples = settling_samples + analysis_samples
    time_axis = np.arange(samples) / SAMPLE_RATE
    audio = np.sin(2.0 * np.pi * frequency * time_axis)
    fold = np.ones(samples)
    symmetry = np.zeros(samples)
    variants = {
        "2x direct": (dsp.wavefolder_external_x2, False),
        "4x direct": (dsp.wavefolder_external_x4, False),
        "4x + ADAA": (dsp.wavefolder_external_x4, True),
        "16x direct": (dsp.wavefolder_external_x16, False),
    }
    results = []
    for character_name, character in CHARACTERS.items():
        for mode, (function, adaa) in variants.items():
            output = function(audio, fold, symmetry, SAMPLE_RATE, adaa, character)[
                settling_samples:
            ]
            results.append(
                (
                    character_name,
                    mode,
                    off_harmonic_energy_db(output, frequency),
                )
            )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=1.0)
    parser.add_argument("--cpu-samples", type=int, default=480_000)
    parser.add_argument("--cpu-repeats", type=int, default=7)
    arguments = parser.parse_args()

    print("Quality at 48 kHz against a 768 kHz render with FIR decimation")
    print("| Character | Scenario | Mode | Magnitude-spectrum error |")
    print("| --- | --- | ---: | ---: |")
    for character, scenario, mode, error in quality_results(arguments.seconds):
        print(f"| {character} | {scenario} | {mode} | {error:.2f} dB |")

    print()
    print("Static timbre for a 100 Hz sine through the 4x direct path")
    print("| Character | Fold | RMS | Harmonic centroid | H3/H1 | H5/H1 | H7/H1 |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for character, fold, rms, centroid, h3, h5, h7 in timbre_results():
        print(
            f"| {character} | {fold:.2f} | {rms:.3f} | {centroid:.0f} Hz | "
            f"{h3:.1f} dB | {h5:.1f} dB | {h7:.1f} dB |"
        )

    print()
    print("Alias stress: 4.187 kHz external sine, full fold, 48 kHz host")
    print("| Character | Mode | Non-harmonic energy |")
    print("| --- | --- | ---: |")
    for character, mode, alias_energy in alias_results():
        print(f"| {character} | {mode} | {alias_energy:.1f} dB |")

    print()
    print("CPU cost (median Python binding render time)")
    print("| Character | Mode | Time per host sample | Relative to 2x |")
    print("| --- | --- | ---: | ---: |")
    scenario = SCENARIOS[-1]
    for character_name, character in CHARACTERS.items():
        timings = {
            name: cpu_result(
                function,
                adaa,
                character,
                scenario,
                arguments.cpu_samples,
                arguments.cpu_repeats,
            )
            for name, (function, adaa) in VARIANTS.items()
        }
        baseline = timings["2x direct"]
        for name, timing in timings.items():
            print(
                f"| {character_name} | {name} | {timing:.1f} ns | "
                f"{timing / baseline:.2f}x |"
            )


if __name__ == "__main__":
    main()
