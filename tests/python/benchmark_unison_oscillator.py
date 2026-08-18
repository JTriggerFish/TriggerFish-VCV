"""Measure the standalone unison oscillator's quality and CPU scaling."""

from __future__ import annotations

import argparse
import time

import numpy as np

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000


def render(voices: int, samples: int, pulse_mix: float = 0.0):
    frequency = np.full(samples, 440.0)
    phase = np.arange(samples) / SAMPLE_RATE
    pulse_width = 0.5 + 0.4 * np.sin(2.0 * np.pi * 7.3 * phase)
    return dsp.stacked_oscillator(
        frequency,
        pulse_width,
        voices=voices,
        spread_cents=4.0,
        pulse_mix=pulse_mix,
        width=0.65,
        sample_rate=SAMPLE_RATE,
    )


def timing(voices: int, samples: int, repeats: int, pulse_mix: float):
    render(voices, 1024, pulse_mix)
    values = []
    checksum = 0.0
    for _ in range(repeats):
        started = time.perf_counter_ns()
        output = render(voices, samples, pulse_mix)
        values.append(time.perf_counter_ns() - started)
        checksum += float(output[-1, 0])
    if not np.isfinite(checksum):
        raise RuntimeError("non-finite benchmark output")
    per_sample = float(np.median(values)) / samples
    return per_sample, per_sample * SAMPLE_RATE / 1e9 * 100.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=48_000)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()

    print("voices  waveform  ns/host-sample  one-channel CPU  16ch estimate")
    for voices in (1, 4, 7, 16):
        for waveform, mix in (("saw", 0.0), ("PWM", 1.0)):
            nanoseconds, cpu = timing(voices, args.samples, args.repeats, mix)
            print(
                f"{voices:>6}  {waveform:>8}  {nanoseconds:>14.1f}  "
                f"{cpu:>14.3f}%  {16.0 * cpu:>12.3f}%"
            )


if __name__ == "__main__":
    main()
