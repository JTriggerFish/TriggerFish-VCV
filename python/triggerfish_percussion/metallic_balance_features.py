"""Causal band filters, centered power smoothing and resolved spectra."""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, sosfilt

from .power_envelope import smoothed_power
from .transforms import StftConfig, stft

ATTACK_BINS = ((0, 0.001), (0.001, 0.003), (0.003, 0.01), (0.01, 0.03), (0.03, 0.1))


class MetallicBalanceFeatures:
    def __init__(self, rate, frames):
        if not np.isfinite(rate) or rate < 1000 or frames < 1:
            raise ValueError("Expected a positive audio sample rate and frame count")
        self.rate, self.frames = rate, frames
        duration = frames / rate
        if duration < 3.1:
            raise ValueError(
                "Metallic trajectory fitting requires at least 3.1 seconds"
            )
        self.regions = ((0, 0.12), (0.12, 0.5), (0.5, 1.5), (1.5, 3), (3, duration))
        self.hop = max(1, round(0.008 * rate))
        upper = min(16000, 0.49 * rate)
        edges = tuple(
            hz for hz in (40, 125, 300, 700, 1500, 3000, 6000) if hz < upper
        ) + (upper,)
        self.bands = tuple((a, b) for a, b in zip(edges, edges[1:]) if a < b)
        self.filters = [
            butter(2, band, btype="bandpass", fs=rate, output="sos")
            for band in self.bands
        ]
        self.times = np.arange(0, frames, self.hop) / rate

    def __call__(self, samples):
        if np.iscomplexobj(samples):
            raise ValueError("Candidate must be real mono audio")
        samples = np.asarray(samples, dtype=np.float64)
        if samples.shape != (self.frames,) or not np.isfinite(samples).all():
            raise ValueError(
                "Candidate must be finite and match the reference duration"
            )
        envelopes = np.array(
            [
                smoothed_power(sosfilt(sos, samples), self.rate * 0.012)[:: self.hop]
                for sos in self.filters
            ]
        )
        resolved = stft(samples, self.rate, StftConfig(16384, 4096))
        selected = (resolved.frequencies_hz >= 40) & (
            resolved.frequencies_hz <= min(16000, 0.49 * self.rate)
        )
        self.frequencies = resolved.frequencies_hz[selected]
        power = resolved.power[selected]
        spectra = []
        for a, b in self.regions:
            selected_times = (resolved.times_seconds >= a) & (
                resolved.times_seconds < b
            )
            if not selected_times.any():
                raise ValueError(
                    f"Resolved STFT has no frames in region {a:g}–{b:g} seconds; "
                    "sample rate/duration do not support this measurement"
                )
            spectra.append(power[:, selected_times].mean(axis=1))
        spectra = np.array(spectra)
        # Sub-bin smoothing only: keep narrow ridges, average stochastic beating
        # over the region rather than smearing across adjacent resonances.
        spectra = gaussian_filter1d(spectra, 0.5, axis=1)
        attack = stft(
            samples[: round(0.14 * self.rate)], self.rate, StftConfig(512, 128)
        )
        selected = (attack.frequencies_hz >= 80) & (
            attack.frequencies_hz <= min(16000, 0.49 * self.rate)
        )
        result = dict(
            envelope=envelopes,
            spectrum=spectra,
            attack=attack.power[selected][:, attack.times_seconds < 0.12],
            transient=np.array(
                [
                    np.mean(samples[round(a * self.rate) : round(b * self.rate)] ** 2)
                    for a, b in ATTACK_BINS
                ]
            ),
        )
        if any(not np.isfinite(value).all() for value in result.values()):
            raise ValueError("Nonfinite metallic measurement")
        return result
