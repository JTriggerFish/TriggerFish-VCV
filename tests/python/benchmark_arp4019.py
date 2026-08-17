"""Print reproducible ARP 4019 model/reference comparisons."""

from __future__ import annotations

import numpy as np

import _triggerfish_dsp as dsp
from reference_arp4019 import circuit_values

SAMPLE_RATE = 48_000.0


def harmonic_amplitudes(signal, frequency, sample_rate, count=9, cycles=100):
    selected_count = int(round(cycles * sample_rate / frequency))
    selected = np.asarray(signal[-selected_count:])
    time = np.arange(selected_count) / sample_rate
    return np.asarray(
        [
            2.0
            * abs(np.mean(selected * np.exp(-2j * np.pi * harmonic * frequency * time)))
            for harmonic in range(1, count + 1)
        ]
    )


def continuous_proxy(amplitude, frequency, duration=0.25, factor=16):
    values = circuit_values()
    sample_rate = factor * SAMPLE_RATE
    time = np.arange(int(duration * sample_rate)) / sample_rate
    audio = amplitude * np.sin(2.0 * np.pi * frequency * time)
    differential = values["audio_input_scale"] * audio
    current = values["unity_control_current_amps"] * np.tanh(
        differential / (2.0 * 25.85e-3)
    )
    target = values["output_feedback_resistance_ohms"] * current
    coefficient = -np.expm1(-2.0 * np.pi * values["output_bandwidth_hz"] / sample_rate)
    output = np.empty_like(target)
    state = 0.0
    for index, value in enumerate(target):
        state += coefficient * (value - state)
        output[index] = state
    return output, sample_rate


def main():
    values = circuit_values()
    print("small-signal x4 magnitude error against the analog 56k/100p pole")
    print("frequency Hz | error dB")
    for frequency in (1000.0, 5000.0, 10000.0, 15000.0, 20000.0):
        count = int(0.2 * SAMPLE_RATE)
        time = np.arange(count) / SAMPLE_RATE
        audio = 1.0e-3 * np.sin(2.0 * np.pi * frequency * time)
        output = dsp.arp4019_x4(audio, np.full(count, 10.0), np.full(count, -10.0))
        measured = harmonic_amplitudes(output, frequency, SAMPLE_RATE, 1)[0] / 1e-3
        expected = 1.0 / np.sqrt(1.0 + (frequency / values["output_bandwidth_hz"]) ** 2)
        print(f"{frequency:12.0f} | {20 * np.log10(measured / expected):+8.4f}")

    print("\n1 kHz nonlinear transfer against a 16x continuous-time proxy")
    print("input Vpk | model THD | reference THD")
    for amplitude in (0.1, 1.0, 5.0, 10.0):
        count = int(0.25 * SAMPLE_RATE)
        time = np.arange(count) / SAMPLE_RATE
        audio = amplitude * np.sin(2.0 * np.pi * 1000.0 * time)
        model = dsp.arp4019_x4(audio, np.full(count, 10.0), np.full(count, -10.0))
        reference, reference_rate = continuous_proxy(amplitude, 1000.0)
        model_h = harmonic_amplitudes(model, 1000.0, SAMPLE_RATE)
        reference_h = harmonic_amplitudes(reference, 1000.0, reference_rate)
        model_thd = np.linalg.norm(model_h[1:]) / model_h[0]
        reference_thd = np.linalg.norm(reference_h[1:]) / reference_h[0]
        print(
            f"{amplitude:9.1f} | {100 * model_thd:9.6f}% | "
            f"{100 * reference_thd:12.6f}%"
        )


if __name__ == "__main__":
    main()
