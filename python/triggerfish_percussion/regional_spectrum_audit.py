"""Inspect missing/excess spectral energy without summing it into a fit score.

Equal log-frequency cells cover the entire audible range, not protected notes.
Welch windows stay inside each region: late bloom cannot leak into the attack.
The reference fixes both the analysis floor and the missing-energy mask.
"""

import numpy as np
from scipy.signal import welch


class RegionalSpectrumAudit:
    def __init__(self, reference, rate, regions, bins_per_octave=12):
        reference = np.asarray(reference)
        self.rate, self.frames = rate, reference.size
        self.regions = np.asarray(regions, dtype=float)
        if (
            reference.ndim != 1
            or not np.isfinite(reference).all()
            or not np.isfinite(rate)
            or rate < 32000
            or self.regions.ndim != 2
            or self.regions.shape[1] != 2
            or not len(self.regions)
            or not np.isfinite(self.regions).all()
            or np.any(self.regions[:, 0] < 0)
            or np.any(self.regions[:, 1] > self.frames / rate)
            or np.any(np.diff(self.regions, axis=1) < 0.05 - 1e-12)
            or not 1 <= bins_per_octave <= 48
        ):
            raise ValueError("Expected finite mono audio and valid >=50 ms regions")
        # A uniform rule throughout the spectrum; no gong-specific boundaries.
        count = int(np.ceil(np.log2(15000 / 20) * bins_per_octave))
        self.edges = np.geomspace(20, 15000, count + 1)
        target = self.power(reference)
        self.floor = np.maximum(target.max(axis=1, keepdims=True), 1e-20) * 1e-6
        self.target = self.db(target)
        self.active = target > self.floor * 1000  # Reference regional peak -30 dB.

    def power(self, samples):
        samples = np.asarray(samples)
        if samples.shape != (self.frames,) or not np.isfinite(samples).all():
            raise ValueError("Expected finite matching mono audio")
        result = []
        for start, end in self.regions:
            segment = samples[round(start * self.rate) : round(end * self.rate)]
            f, p = welch(
                segment, self.rate, nperseg=min(8192, len(segment)), nfft=32768
            )
            # Integrate via the cumulative PSD; narrow cells never become empty.
            df = f[1] - f[0]
            integral = np.r_[0, np.cumsum(p) * df]
            boundaries = np.r_[f - df / 2, f[-1] + df / 2]
            result.append(np.diff(np.interp(self.edges, boundaries, integral)))
        return np.array(result)

    def db(self, power):
        return 10 * np.log10(np.maximum(power, self.floor))

    def measure(self, samples):
        power = self.power(samples)
        error = self.db(power) - self.target
        active = self.active | (power > self.floor * 1000)
        return dict(
            edges_hz=self.edges.tolist(),
            regions_seconds=self.regions.tolist(),
            reference_db=self.target.tolist(),
            error_db=error.tolist(),
            reference_active=self.active.tolist(),
            max_deficit_db=float(
                np.max(np.where(self.active, np.maximum(-error, 0), 0))
            ),
            max_excess_db=float(np.max(np.where(active, np.maximum(error, 0), 0))),
            note="Diagnostic, not a perceptual score or automatic acceptance gate",
        )
