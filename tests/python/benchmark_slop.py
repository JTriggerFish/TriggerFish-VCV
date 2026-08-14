import statistics
import time

import numpy as np

import _triggerfish_dsp as dsp

SIZE = 2_000_000
ROUNDS = 9


def benchmark_group(processors, pitch, detune):
    timings = {name: [] for name in processors}
    names = list(processors)
    for processor in processors.values():
        processor(pitch, detune)

    for round_index in range(ROUNDS):
        offset = round_index % len(names)
        for name in names[offset:] + names[:offset]:
            start = time.perf_counter()
            processors[name](pitch, detune)
            timings[name].append(time.perf_counter() - start)

    return {
        name: SIZE / statistics.median(samples) / 1_000_000.0
        for name, samples in timings.items()
    }


def main():
    phase = np.linspace(0.0, 2.0 * np.pi, SIZE)
    pitch = 5.0 * np.sin(phase)
    detune = 2.0 * np.sin(17.0 * phase)
    processors = {
        "legacy double": dsp.detune_legacy_double,
        "optimized mixed": dsp.detune_optimized,
    }
    results = benchmark_group(processors, pitch, detune)
    baseline = results["legacy double"]

    print("Slop linear-Hz detune conversion")
    for name, throughput in results.items():
        change = 100.0 * (throughput / baseline - 1.0)
        four_channel_streams = throughput * 1_000_000.0 / (4.0 * 48_000.0)
        print(
            f"  {name:18s} {throughput:7.2f} Mconversions/s "
            f"{change:+6.2f}%  ({four_channel_streams:6.1f} equivalent "
            "four-channel conversion streams/core at 48 kHz)"
        )


if __name__ == "__main__":
    main()
