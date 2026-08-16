"""Compare oscillator quality and CPU cost at the available sample rates."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import time

import numpy as np
from scipy.signal import resample_poly

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000
REFERENCE_FACTOR = 16
REFERENCE_WINDOW = ("kaiser", 12.0)
PRODUCTION_VARIANTS = {
    "2x": dsp.tb303_oscillator_x2,
    "4x": dsp.tb303_oscillator_x4,
}


@dataclass(frozen=True)
class Scenario:
    name: str
    frequency: float
    wave: float
    fm_amplitude: float = 0.0
    fm_frequency: float = 0.0
    linear_fm: bool = False
    sync_frequency: float | None = None


SCENARIOS = (
    Scenario("1 kHz square", 1_000.0, 1.0),
    Scenario("6 kHz saw", 6_000.0, 0.0),
    Scenario("6 kHz square", 6_000.0, 1.0),
    Scenario("12 kHz square", 12_000.0, 1.0),
    Scenario("6.1 kHz saw, 997 Hz sync", 6_100.0, 0.0, sync_frequency=997.0),
    Scenario("1.5 kHz saw, exponential FM", 1_500.0, 0.0, 4.0, 997.0),
    Scenario("1.2 kHz saw, through-zero FM", 1_200.0, 0.0, 8.0, 733.0, True),
    Scenario(
        "1.2 kHz morph, through-zero FM + sync",
        1_200.0,
        0.5,
        8.0,
        733.0,
        True,
        997.0,
    ),
)


def inputs(scenario: Scenario, sample_rate: int, seconds: float):
    samples = round(seconds * sample_rate)
    time_axis = np.arange(samples, dtype=np.float64) / sample_rate
    pitch = np.full(samples, np.log2(scenario.frequency / 261.625565))
    zero = np.zeros(samples)
    fm = scenario.fm_amplitude * np.sin(
        2.0 * np.pi * scenario.fm_frequency * time_axis + 0.19
    )
    wave = np.full(samples, scenario.wave)
    sync = None
    if scenario.sync_frequency is not None:
        sync = 5.0 * np.sin(2.0 * np.pi * scenario.sync_frequency * time_axis + 0.37)
    return pitch, zero, fm, wave, sync


def render(function, scenario: Scenario, sample_rate: int, seconds: float):
    pitch, zero, fm, wave, sync = inputs(scenario, sample_rate, seconds)
    return function(
        pitch,
        zero,
        fm,
        zero,
        wave,
        sample_rate=sample_rate,
        linear_fm=scenario.linear_fm,
        sync=sync,
    )[:, 2]


def aligned_error(candidate, reference, maximum_lag=128):
    """Gain- and integer-delay-aligned waveform residual in dB."""
    candidate = candidate - candidate.mean()
    reference = reference - reference.mean()
    best = np.inf
    for lag in range(-maximum_lag, maximum_lag + 1):
        if lag >= 0:
            left = candidate[lag:]
            right = reference[: left.size]
        else:
            right = reference[-lag:]
            left = candidate[: right.size]
        gain = np.dot(left, right) / np.dot(left, left)
        error = np.sqrt(np.mean(np.square(gain * left - right)))
        best = min(best, error / np.sqrt(np.mean(np.square(right))))
    return 20.0 * np.log10(best)


def spectral_error(candidate, reference):
    """Normalized magnitude-spectrum error in dB; more negative is better."""
    window = np.hanning(reference.size)
    reference_spectrum = np.abs(np.fft.rfft(reference * window))
    candidate_spectrum = np.abs(np.fft.rfft(candidate * window))
    reference_spectrum /= np.linalg.norm(reference_spectrum)
    candidate_spectrum /= np.linalg.norm(candidate_spectrum)
    return 20.0 * np.log10(np.linalg.norm(candidate_spectrum - reference_spectrum))


def quality_results(seconds: float):
    results = []
    high_rate = SAMPLE_RATE * REFERENCE_FACTOR
    padding = 0.25
    render_seconds = seconds + 2.0 * padding
    for scenario in SCENARIOS:
        reference = resample_poly(
            render(dsp.tb303_oscillator_x1, scenario, high_rate, render_seconds),
            1,
            REFERENCE_FACTOR,
            window=REFERENCE_WINDOW,
        )
        start = round(padding * SAMPLE_RATE)
        stop = start + round(seconds * SAMPLE_RATE)
        reference = reference[start:stop]
        for variant, function in PRODUCTION_VARIANTS.items():
            candidate = render(function, scenario, SAMPLE_RATE, render_seconds)
            candidate = candidate[start:stop]
            results.append(
                (
                    scenario.name,
                    variant,
                    spectral_error(candidate, reference),
                    aligned_error(candidate, reference),
                )
            )
    return results


def cpu_result(function, scenario: Scenario, samples: int, repeats: int):
    seconds = samples / SAMPLE_RATE
    arguments = inputs(scenario, SAMPLE_RATE, seconds)
    pitch, zero, fm, wave, sync = arguments

    def run():
        return function(
            pitch,
            zero,
            fm,
            zero,
            wave,
            sample_rate=SAMPLE_RATE,
            linear_fm=scenario.linear_fm,
            sync=sync,
        )

    run()
    timings = []
    checksum = 0.0
    for _ in range(repeats):
        started = time.perf_counter_ns()
        output = run()
        timings.append(time.perf_counter_ns() - started)
        checksum += float(output[-1, 2])
    if not np.isfinite(checksum):
        raise RuntimeError("non-finite benchmark output")
    return float(np.median(timings)) / samples


def print_quality_table(results):
    print("Quality at 48 kHz against a 768 kHz render with FIR decimation")
    print("Spectral error includes aliasing and passband-magnitude differences.")
    print("| Scenario | Mode | Spectral error | Aligned waveform residual |")
    print("| --- | ---: | ---: | ---: |")
    for scenario, variant, spectrum, waveform in results:
        print(f"| {scenario} | {variant} | {spectrum:.2f} dB | {waveform:.2f} dB |")


def print_cpu_table(samples: int, repeats: int):
    cases = (SCENARIOS[0], SCENARIOS[-1])
    print()
    print(f"CPU cost, median of {repeats} renders ({samples:,} samples each)")
    print("| Scenario | Mode | Time per sample | Relative to 2x |")
    print("| --- | ---: | ---: | ---: |")
    for scenario in cases:
        timings = {
            name: cpu_result(function, scenario, samples, repeats)
            for name, function in PRODUCTION_VARIANTS.items()
        }
        for name, timing in timings.items():
            print(
                f"| {scenario.name} | {name} | {timing:.1f} ns | "
                f"{timing / timings['2x']:.2f}x |"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Keep the default invocation identical to the published technical-report
    # run. Shorter exploratory runs can still override each value explicitly.
    parser.add_argument("--seconds", type=float, default=2.0)
    parser.add_argument("--cpu-samples", type=int, default=960_000)
    parser.add_argument("--cpu-repeats", type=int, default=9)
    arguments = parser.parse_args()
    print_quality_table(quality_results(arguments.seconds))
    print_cpu_table(arguments.cpu_samples, arguments.cpu_repeats)


if __name__ == "__main__":
    main()
