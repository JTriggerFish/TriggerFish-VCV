"""Controlled faults, applied to a reference without candidate normalization."""

import numpy as np
from scipy.signal import butter, sosfilt


def fault_examples(reference, rate):
    """Return three severities per fault; these are tests, not decompositions.

    Pitch is shifted by resampling (and therefore changes duration too); label
    that confound explicitly instead of calling it an isolated pitch test.
    """
    from scipy.signal import resample_poly

    time = np.arange(len(reference)) / rate
    attack_rms = np.sqrt(np.mean(reference[: round(0.03 * rate)] ** 2))
    # A zero-phase shelf isolates attenuation. Causal low-pass subtraction
    # can instead BOOST bins through phase cancellation.
    frequencies = np.fft.rfftfreq(len(reference), 1 / rate)
    bass_weight = 1 / (1 + (frequencies / 180) ** 8)
    bass = np.fft.irfft(np.fft.rfft(reference) * bass_weight, n=len(reference))
    noise = sosfilt(
        butter(4, [700, 6000], fs=rate, btype="bandpass", output="sos"),
        np.random.default_rng(735).standard_normal(len(reference)),
    )
    noise /= np.sqrt(np.mean(noise[: round(0.03 * rate)] ** 2))
    ring = np.sin(2 * np.pi * 700 * time) * np.exp(-np.log(1000) * time / 0.4)
    for severity in (1, 2, 3):
        yield "missing-bass", severity, reference - bass * (0.2 * severity)
        yield "excess-ringing", severity, reference + ring * attack_rms * (
            0.2 * severity
        )
        envelope = np.exp(-np.log(1000) * time / (0.08 * severity))
        yield "long-noise", severity, reference + 0.4 * attack_rms * noise * envelope
        ratio = 2 ** (severity / 12)
        shifted = resample_poly(reference, 1000, round(1000 * ratio))
        yield "pitch-and-speed-up", severity, np.pad(
            shifted, (0, len(reference) - len(shifted))
        )
