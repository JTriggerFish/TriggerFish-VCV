"""Fixed-level STFT-band envelopes with explicit onset-to-bloom contrast.

Unlike causal low-order band filters, frequency masks do not let a strong low
partial masquerade as a delayed high-band signal. Windows smear time equally
on reference and candidate; this is an offline analysis, not a causal filter.
"""

import numpy as np
from scipy.signal import stft


class SpectralBloomLoss:
    units = "dB spectral-envelope and rise error"
    influence_threshold = 0.01

    def __init__(self, reference, rate):
        self.rate, self.frames = rate, len(reference)
        self.edges = np.geomspace(80, min(16000, rate * 0.45), 25)
        self.regions = [
            (a, b)
            for a, b in zip(
                (0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1, 1.3, 1.7, 2.2, 3, 4, 5),
                (0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1, 1.3, 1.7, 2.2, 3, 4, 5, 6),
            )
        ]
        target = self.power(reference)
        self.floor = max(target.max() * 1e-7, 1e-20)
        self.target = self.db(target)
        # A whole-sound peak would suppress the quiet pre-bloom cells. Keep
        # those cells and ask explicitly for their contrast with the rise.
        self.active = target.max(axis=1) > target.max() * 1e-5
        self.specification = dict(
            version="spectral-bloom-v2",
            fft_size=4096,
            hop_seconds=0.01,
            bands_hz=self.edges.tolist(),
            regions=self.regions,
            floor_below_reference_peak_db=70,
            normalization=False,
            rise_contrast_weight=1,
            spectral_envelope_weight=1,
        )

    def power(self, samples):
        samples = np.asarray(samples)
        if samples.shape != (self.frames,) or not np.isfinite(samples).all():
            raise ValueError("Expected finite mono audio matching the reference")
        f, t, z = stft(
            samples,
            self.rate,
            nperseg=4096,
            noverlap=4096 - round(0.01 * self.rate),
            boundary="zeros",
        )
        bands = np.array(
            [
                np.sum(abs(z[(f >= lo) & (f < hi)]) ** 2, axis=0)
                for lo, hi in zip(self.edges[:-1], self.edges[1:])
            ]
        )
        return np.array(
            [bands[:, (t >= lo) & (t < hi)].mean(axis=1) for lo, hi in self.regions]
        ).T

    def db(self, power):
        return 10 * np.log10(np.maximum(power, self.floor))

    def residual(self, samples, regions=None):
        error = (self.db(self.power(samples)) - self.target)[self.active]
        rise = error[:, 2:9] - error[:, :1]
        return np.concatenate(
            (error.ravel() / np.sqrt(error.size), rise.ravel() / np.sqrt(rise.size))
        )

    def diagnostics(self, samples):
        error = (self.db(self.power(samples)) - self.target)[self.active]
        return dict(
            envelope_rms_db=float(np.sqrt(np.mean(error**2))),
            rise_rms_db=float(np.sqrt(np.mean((error[:, 2:9] - error[:, :1]) ** 2))),
            difference_db=error.tolist(),
        )
