"""Compare oscillator oversampling filters with a 16x/FIR reference."""

from __future__ import annotations

import numpy as np
from scipy.signal import resample_poly

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000
SECONDS = 1.0
VARIANTS = {
    "1x": dsp.tb303_oscillator_x1,
    "2x Butterworth-5": dsp.tb303_oscillator_x2_order5,
    "2x Chebyshev-7": dsp.tb303_oscillator_x2,
    "4x Butterworth-5": dsp.tb303_oscillator_x4_order5,
    "4x Chebyshev-7": dsp.tb303_oscillator_x4,
}


def render(function, frequency, sample_rate):
    samples = round(SECONDS * sample_rate)
    pitch = np.full(samples, np.log2(frequency / 261.625565))
    zero = np.zeros(samples)
    square = np.ones(samples)
    return function(pitch, zero, zero, zero, square, sample_rate=sample_rate)[:, 2]


def aligned_error(candidate, reference, maximum_lag=128):
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


def magnitude_error(candidate, reference):
    window = np.hanning(reference.size)
    reference_spectrum = np.abs(np.fft.rfft(reference * window))
    candidate_spectrum = np.abs(np.fft.rfft(candidate * window))
    reference_spectrum /= np.linalg.norm(reference_spectrum)
    candidate_spectrum /= np.linalg.norm(candidate_spectrum)
    return 20.0 * np.log10(np.linalg.norm(candidate_spectrum - reference_spectrum))


def main():
    for frequency in (1_000.0, 6_000.0, 12_000.0):
        high_rate = SAMPLE_RATE * 16
        reference = resample_poly(
            render(dsp.tb303_oscillator_x1, frequency, high_rate),
            1,
            16,
            window=("kaiser", 12.0),
        )
        reference = reference[reference.size // 2 :]
        print(f"{frequency:g} Hz square")
        for name, function in VARIANTS.items():
            candidate = render(function, frequency, SAMPLE_RATE)
            candidate = candidate[candidate.size // 2 :]
            print(
                f"  {name:20s} magnitude {magnitude_error(candidate, reference):7.2f} dB"
                f"  phase-sensitive {aligned_error(candidate, reference):7.2f} dB"
            )


if __name__ == "__main__":
    main()
