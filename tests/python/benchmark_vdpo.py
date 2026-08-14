import statistics
import time

import numpy as np

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000
DURATION = 2.0
ROUNDS = 9
C4 = 261.625565


def benchmark_group(processors, third_input, damping):
    size = int(SAMPLE_RATE * DURATION)
    audio = np.zeros(size)
    damping_values = np.full(size, damping)
    values = np.full(size, third_input)
    timings = {name: [] for name in processors}
    names = list(processors)

    for processor in processors.values():
        processor(audio, damping_values, values, SAMPLE_RATE)

    for round_index in range(ROUNDS):
        offset = round_index % len(names)
        for name in names[offset:] + names[:offset]:
            start = time.perf_counter()
            processors[name](audio, damping_values, values, SAMPLE_RATE)
            timings[name].append(time.perf_counter() - start)

    return {
        name: size / statistics.median(samples) / 1_000_000.0
        for name, samples in timings.items()
    }


def print_results(title, results):
    baseline = next(iter(results.values()))
    print(title)
    for name, throughput in results.items():
        change = 100.0 * (throughput / baseline - 1.0)
        print(f"  {name:20s} {throughput:7.3f} MSamples/s  {change:+6.2f}%")


def main():
    integrators = {
        "legacy BDF": dsp.vdpo_bdf,
        "split integrator": dsp.vdpo,
    }
    pitch_processors = {
        "std::exp2 baseline": dsp.vdpo_pitch_std,
        "Rack fast exp2": dsp.vdpo_pitch_fast_exp2,
    }

    for damping, frequency in ((0.5, 1_000.0), (5.0, 8_000.0), (9.0, 18_000.0)):
        print(f"\ndamping={damping:g}, frequency={frequency:g} Hz")
        print_results(
            "Integrator comparison:",
            benchmark_group(integrators, 2.0 * np.pi * frequency, damping),
        )
        print_results(
            "End-to-end pitch conversion variants:",
            benchmark_group(
                pitch_processors,
                np.log2(frequency / C4),
                damping,
            ),
        )


if __name__ == "__main__":
    main()
