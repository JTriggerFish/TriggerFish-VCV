"""Quality and throughput comparison for diode-ladder resampler candidates."""

from time import perf_counter

import numpy as np
from scipy.signal import resample_poly

import _triggerfish_dsp as dsp

PROCESSORS = {
    "2x order 7": dsp.diode_ladder_x2,
    "2x order 9": dsp.diode_ladder_x2_order9,
    "4x order 7": dsp.diode_ladder_x4,
    "4x order 9": dsp.diode_ladder_x4_order9,
}

ROUND_TRIPS = {
    "2x order 7": dsp.resampler_round_trip_x2_order7,
    "2x order 9": dsp.resampler_round_trip_x2_order9,
    "4x order 7": dsp.resampler_round_trip_x4_order7,
    "4x order 9": dsp.resampler_round_trip_x4_order9,
}


def benchmark(processor, label, signal, repeats=5):
    timings = []
    for _ in range(repeats):
        started = perf_counter()
        processor(signal, 1_200.0, 0.65, False, 1.0, 0.5, 48_000.0)
        timings.append(perf_counter() - started)
    elapsed = min(timings)
    print(
        f"{label:>10}: {signal.size / elapsed / 1e6:7.3f} Msamples/s "
        f"({elapsed / signal.size * 1e9:7.1f} ns/host sample)"
    )


def oversampling_comparison(frequency, cutoff, resonance, drive):
    sample_rate = 48_000.0
    size = int(2.0 * sample_rate)
    time = np.arange(size) / sample_rate
    audio = 5.0 * np.sin(2.0 * np.pi * frequency * time)
    print(f"\n{frequency:.0f} Hz input, {cutoff:.0f} Hz cutoff, {drive:.1f}x drive")
    for order, x2_processor, x4_processor in (
        (7, dsp.diode_ladder_x2, dsp.diode_ladder_x4),
        (9, dsp.diode_ladder_x2_order9, dsp.diode_ladder_x4_order9),
    ):
        x2 = x2_processor(audio, cutoff, resonance, False, drive, 0.5, sample_rate)[
            size // 2 :
        ]
        x4 = x4_processor(audio, cutoff, resonance, False, drive, 0.5, sample_rate)[
            size // 2 :
        ]
        # Align the candidates by cross-correlation before comparing waveforms.
        correlation = np.correlate(x4[:4096], x2[:4096], mode="full")
        lag = int(np.argmax(correlation) - 4095)
        if lag > 0:
            x4, x2 = x4[lag:], x2[:-lag]
        elif lag < 0:
            x4, x2 = x4[:lag], x2[-lag:]
        relative_error = np.sqrt(np.mean((x2 - x4) ** 2)) / np.sqrt(np.mean(x4**2))
        print(
            f"  order {order}: 2x/4x difference {100.0 * relative_error:6.2f}% RMS "
            f"(integer lag {lag:+d})"
        )


def alias_residual(processor, frequency):
    """Measure non-harmonic output from a coherent periodic sine.

    A time-invariant analog nonlinearity driven by a settled sine can only produce
    integer harmonics. Energy at other host-rate FFT bins is therefore numerical
    noise, modulation, or folded energy from harmonics above Nyquist.
    """

    sample_rate = 48_000.0
    size = int(4.0 * sample_rate)
    time = np.arange(size) / sample_rate
    audio = 5.0 * np.sin(2.0 * np.pi * frequency * time)
    output = processor(audio, 20_000.0, 0.0, False, 66.6, 0.0, sample_rate)[
        -int(sample_rate) :
    ].astype(np.float64)
    spectrum = np.fft.rfft(output)
    harmonic_mask = np.zeros(spectrum.size, dtype=bool)
    harmonic_mask[0] = True
    for harmonic in range(1, int(sample_rate // (2.0 * frequency)) + 1):
        bin_index = int(harmonic * frequency)
        harmonic_mask[max(0, bin_index - 1) : bin_index + 2] = True
    wanted = np.linalg.norm(spectrum[harmonic_mask])
    residual = np.linalg.norm(spectrum[~harmonic_mask])
    ratio_db = 20.0 * np.log10(max(residual / wanted, np.finfo(float).tiny))
    return ratio_db


def high_rate_alias_reference(frequency):
    """Render at 16x and FIR-decimate to estimate the bandlimited reference."""

    sample_rate = 48_000.0
    factor = 16
    reference_rate = factor * sample_rate
    size = int(4.0 * reference_rate)
    time = np.arange(size) / reference_rate
    audio = 5.0 * np.sin(2.0 * np.pi * frequency * time)
    output = dsp.diode_ladder_x1(audio, 20_000.0, 0.0, False, 66.6, 0.0, reference_rate)
    downsampled = resample_poly(output, 1, factor, window=("kaiser", 14.0))
    spectrum = np.fft.rfft(downsampled[-int(sample_rate) :].astype(np.float64))
    harmonic_mask = np.zeros(spectrum.size, dtype=bool)
    harmonic_mask[0] = True
    for harmonic in range(1, int(sample_rate // (2.0 * frequency)) + 1):
        bin_index = int(harmonic * frequency)
        harmonic_mask[max(0, bin_index - 1) : bin_index + 2] = True
    wanted = np.linalg.norm(spectrum[harmonic_mask])
    residual = np.linalg.norm(spectrum[~harmonic_mask])
    return 20.0 * np.log10(max(residual / wanted, np.finfo(float).tiny))


def vca_alias_residual(frequency, oversampled):
    sample_rate = 48_000.0
    size = int(4.0 * sample_rate)
    time = np.arange(size) / sample_rate
    audio = 5.0 * np.sin(2.0 * np.pi * frequency * time)
    output = dsp.diode_ladder_vca_x4(
        audio,
        np.ones(size),
        20_000.0,
        0.0,
        False,
        66.6,
        0.0,
        sample_rate,
        oversampled,
    )[-int(sample_rate) :]
    spectrum = np.fft.rfft(output.astype(np.float64))
    harmonic_mask = np.zeros(spectrum.size, dtype=bool)
    harmonic_mask[0] = True
    for harmonic in range(1, int(sample_rate // (2.0 * frequency)) + 1):
        bin_index = int(harmonic * frequency)
        harmonic_mask[max(0, bin_index - 1) : bin_index + 2] = True
    wanted = np.linalg.norm(spectrum[harmonic_mask])
    residual = np.linalg.norm(spectrum[~harmonic_mask])
    return 20.0 * np.log10(max(residual / wanted, np.finfo(float).tiny))


def round_trip_phase():
    """Report magnitude and group delay of an identity oversampling round trip."""

    sample_rate = 48_000.0
    fft_size = 1 << 18
    impulse = np.zeros(fft_size)
    impulse[0] = 1.0
    probe_frequencies = np.array([100, 1_000, 5_000, 10_000, 15_000, 20_000])
    bins = np.rint(probe_frequencies * fft_size / sample_rate).astype(int)
    angular_frequency = 2.0 * np.pi * np.fft.rfftfreq(fft_size)
    print("\nIdentity upsample/downsample round trip")
    print(
        "             magnitude at 20 kHz    group delay [samples] at 0.1/1/5/10/15/20 kHz"
    )
    for label, processor in ROUND_TRIPS.items():
        response = np.fft.rfft(processor(impulse))
        phase = np.unwrap(np.angle(response))
        group_delay = -np.gradient(phase, angular_frequency)
        magnitude_db = 20.0 * np.log10(abs(response[bins[-1]]))
        delays = "/".join(f"{group_delay[index]:.2f}" for index in bins)
        print(f"{label:>10}: {magnitude_db:8.3f} dB          {delays}")


if __name__ == "__main__":
    rng = np.random.default_rng(303)
    audio = rng.uniform(-5.0, 5.0, 480_000)
    print("Throughput (secondary)")
    for candidate, processor in PROCESSORS.items():
        benchmark(processor, candidate, audio)

    print("\nHigh-drive non-harmonic residual [dBc]")
    frequencies = (5_000.0, 7_000.0, 9_000.0, 11_000.0)
    print(
        "             "
        + "".join(f"{frequency / 1000:>9.0f} kHz" for frequency in frequencies)
    )
    for candidate, processor in PROCESSORS.items():
        residuals = "".join(
            f"{alias_residual(processor, frequency):>13.2f}"
            for frequency in frequencies
        )
        print(f"{candidate:>10}: {residuals}")
    reference_residuals = "".join(
        f"{high_rate_alias_reference(frequency):>13.2f}" for frequency in frequencies
    )
    print(f"{'16x + FIR':>10}: {reference_residuals}")

    print("\nEffect of moving the built-in VCA into the 4x domain [dBc]")
    for oversampled in (False, True):
        residuals = "".join(
            f"{vca_alias_residual(frequency, oversampled):>13.2f}"
            for frequency in frequencies
        )
        label = "VCA at 4x" if oversampled else "host VCA"
        print(f"{label:>10}: {residuals}")

    round_trip_phase()
    oversampling_comparison(997.0, 1_200.0, 0.65, 1.0)
    oversampling_comparison(5_000.0, 6_000.0, 0.8, 20.0)
    oversampling_comparison(9_000.0, 10_000.0, 0.8, 66.6)
