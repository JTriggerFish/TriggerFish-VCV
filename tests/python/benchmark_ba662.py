"""Compare the reduced 303 OTA VCA with the pre-existing VCA cores.

This is a development benchmark, not a claim that either legacy TriggerFish
core is a BA662 reference. It reports level, distortion, residual after gain
matching, and throughput so topology differences are not confused with simple
loudness differences.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.signal import resample_poly

import _triggerfish_dsp as dsp
from reference_ba662 import static_tb303_vca

SAMPLE_RATE = 48_000.0
DURATION = 4.0
FREQUENCY = 997.0


def rms(signal):
    return float(np.sqrt(np.mean(np.square(signal))))


def harmonic_distortion(signal, frequency=FREQUENCY):
    window = np.hanning(signal.size)
    spectrum = np.abs(np.fft.rfft(signal * window))
    frequencies = np.fft.rfftfreq(signal.size, 1.0 / SAMPLE_RATE)
    fundamental = spectrum[np.argmin(np.abs(frequencies - frequency))]
    harmonics = []
    for harmonic in range(2, 10):
        target = harmonic * frequency
        if target >= SAMPLE_RATE / 2.0:
            break
        harmonics.append(spectrum[np.argmin(np.abs(frequencies - target))])
    return float(np.sqrt(np.sum(np.square(harmonics))) / fundamental)


def timed(function, repeats=8):
    best = float("inf")
    for _ in range(repeats):
        start = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - start)
    return best


def main():
    count = int(SAMPLE_RATE * DURATION)
    time_axis = np.arange(count, dtype=np.float64) / SAMPLE_RATE
    audio64 = 5.0 * np.sin(2.0 * np.pi * FREQUENCY * time_axis)
    control64 = np.ones(count, dtype=np.float64)
    accent64 = np.zeros(count, dtype=np.float64)
    audio32 = audio64.astype(np.float32)
    control32 = control64.astype(np.float32)

    candidates = {
        "reduced BA662": lambda: dsp.tb303_vca(
            audio64, control64, accent64, SAMPLE_RATE
        ),
        "legacy OTA": lambda: dsp.vca_ota_legacy(audio32, control32, SAMPLE_RATE, 1.0),
        "shipped transistor": lambda: dsp.vca_transistor(
            audio32, control32, SAMPLE_RATE, 1.0
        ),
    }
    reference = candidates["reduced BA662"]()[int(SAMPLE_RATE) :]
    reference_rms = rms(reference)

    print("candidate             rms V    THD %   matched residual %   Msamples/s")
    for name, render in candidates.items():
        output = np.asarray(render(), dtype=np.float64)[int(SAMPLE_RATE) :]
        output_rms = rms(output)
        gain = reference_rms / output_rms if output_rms > 0.0 else 0.0
        residual = rms(gain * output - reference) / reference_rms
        elapsed = timed(render)
        print(
            f"{name:20s} {output_rms:7.4f} "
            f"{100.0 * harmonic_distortion(output):8.4f} "
            f"{100.0 * residual:20.3f} {count / elapsed / 1.0e6:12.3f}"
        )

    # Compare the 1x memoryless production transfer with an 8x render before
    # downsampling. This isolates nonlinear aliasing at a deliberately high
    # input frequency; the sub-audio output coupling has negligible effect.
    high_frequency = 9000.0
    oversampling = 8
    high_input = 5.0 * np.sin(2.0 * np.pi * high_frequency * time_axis)
    one_x = dsp.tb303_vca(high_input, control64, accent64, SAMPLE_RATE)
    high_time = np.arange(count * oversampling) / (SAMPLE_RATE * oversampling)
    high_input_oversampled = 5.0 * np.sin(2.0 * np.pi * high_frequency * high_time)
    high_reference = resample_poly(
        static_tb303_vca(
            high_input_oversampled,
            np.ones_like(high_input_oversampled),
            np.zeros_like(high_input_oversampled),
        ),
        1,
        oversampling,
    )
    settled = slice(int(SAMPLE_RATE), count)
    match = rms(high_reference[settled]) / rms(one_x[settled])
    alias_residual = rms(match * one_x[settled] - high_reference[settled]) / rms(
        high_reference[settled]
    )
    print(f"\n9 kHz nominal 1x/8x-reference residual: {100 * alias_residual:.3f}%")


if __name__ == "__main__":
    main()
