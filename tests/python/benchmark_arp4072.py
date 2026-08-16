"""Print reproducible ARP 4072 production/reference comparisons."""

from __future__ import annotations

import argparse

import numpy as np

import _triggerfish_dsp as dsp
from reference_arp4072 import (
    harmonic_amplitudes,
    positive_zero_crossing_frequency,
    render_continuous_self_oscillation,
    render_continuous_sine,
    transfer,
)

SAMPLE_RATE = 48_000.0


def sine_gain(processor, frequency, cutoff, resonance, amplitude=1.0e-3):
    duration = 0.5
    time = np.arange(int(SAMPLE_RATE * duration)) / SAMPLE_RATE
    signal = amplitude * np.sin(2.0 * np.pi * frequency * time)
    output = processor(signal, cutoff, resonance)
    selected = time >= duration / 2.0
    carrier = np.exp(-2j * np.pi * frequency * time[selected])
    return 2.0 * np.mean(output[selected] * carrier) / (-1j * amplitude)


def print_small_signal_table():
    print("small-signal magnitude error against the continuous transfer")
    print("cutoff Hz | input Hz | 2x dB | 4x dB")
    for cutoff, frequency in (
        (8000.0, 8000.0),
        (12000.0, 12000.0),
        (8000.0, 18000.0),
        (12000.0, 18000.0),
    ):
        expected = transfer(frequency, cutoff, 0.5)
        errors = []
        for processor in (dsp.arp4072_x2, dsp.arp4072_x4):
            measured = sine_gain(processor, frequency, cutoff, 0.5)
            errors.append(20.0 * np.log10(abs(measured / expected)))
        print(
            f"{cutoff:9.0f} | {frequency:8.0f} | {errors[0]:+6.3f} | {errors[1]:+6.3f}"
        )


def print_nonlinear_case(frequency, amplitude, cutoff, resonance):
    duration = 0.08
    time, reference = render_continuous_sine(
        frequency, amplitude, cutoff, resonance, duration=duration
    )
    signal = amplitude * np.sin(2.0 * np.pi * frequency * time)
    production = dsp.arp4072_x4(signal, cutoff, resonance)
    expected = harmonic_amplitudes(reference, frequency, SAMPLE_RATE, 9, 40)
    actual = harmonic_amplitudes(production, frequency, SAMPLE_RATE, 9, 40)

    print(
        f"\nnonlinear case: {frequency:g} Hz, {amplitude:g} V, "
        f"cutoff {cutoff:g} Hz, resonance {resonance:g}"
    )
    print("harmonic | reference V | production V | relative error")
    for harmonic in (1, 3, 5, 7, 9):
        index = harmonic - 1
        error = actual[index] / expected[index] - 1.0
        print(
            f"{harmonic:8d} | {expected[index]:11.7f} | "
            f"{actual[index]:12.7f} | {100.0 * error:+8.3f}%"
        )


def print_self_oscillation():
    cutoff = 800.0
    resonance = 1.0
    duration = 2.0
    _, reference = render_continuous_self_oscillation(
        cutoff, resonance, duration=duration
    )
    silence = np.zeros(int(SAMPLE_RATE * duration))
    silence[0] = 1.0e-6
    production = dsp.arp4072_x4(silence, cutoff, resonance)
    reference = reference[len(reference) // 2 :]
    production = production[len(production) // 2 :]
    print("\nself-oscillation at 800 Hz requested cutoff")
    print(
        "reference: "
        f"{positive_zero_crossing_frequency(reference, SAMPLE_RATE):.3f} Hz, "
        f"{np.max(np.abs(reference)):.5f} V peak"
    )
    print(
        "production: "
        f"{positive_zero_crossing_frequency(production, SAMPLE_RATE):.3f} Hz, "
        f"{np.max(np.abs(production)):.5f} V peak"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the slower two-second self-oscillation reference",
    )
    arguments = parser.parse_args()
    print_small_signal_table()
    print_nonlinear_case(1000.0, 50.0, 2000.0, 0.55)
    print_nonlinear_case(1000.0, 80.0, 8000.0, 0.55)
    if arguments.full:
        print_self_oscillation()


if __name__ == "__main__":
    main()
