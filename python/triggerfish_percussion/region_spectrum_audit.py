"""Log-spaced spectral shape checks complement integrated band energy."""

import numpy as np
from .transforms import StftConfig, stft


class RegionSpectrumAudit:
    regions = ((0, 0.03), (0.03, 0.1), (0.1, 0.25))

    def __init__(self, reference, rate):
        self.rate = rate
        self.edges = np.geomspace(125, min(8000, 0.45 * rate), 25)
        self.target = self.power(reference)
        self.floor = max(float(self.target.max()), 1e-20) * 1e-7
        self.mask = self.target > self.floor * 100

    def power(self, samples):
        size = max(2048, 2 ** int(np.ceil(np.log2(self.rate / 24))))
        value = stft(samples, self.rate, StftConfig(size, size // 8))
        power = value.power
        result = []
        for low, high in zip(self.edges[:-1], self.edges[1:]):
            bins = (value.frequencies_hz >= low) & (value.frequencies_hz < high)
            result.append(
                [
                    power[bins][
                        :, (value.times_seconds >= a) & (value.times_seconds < b)
                    ].mean()
                    for a, b in self.regions
                ]
            )
        return np.array(result)

    def error(self, samples):
        return 10 * np.log10(
            np.maximum(self.power(samples), self.floor)
            / np.maximum(self.target, self.floor)
        )

    def measure(self, samples):
        power = self.power(samples)
        error = 10 * np.log10(
            np.maximum(power, self.floor) / np.maximum(self.target, self.floor)
        )
        evaluated = self.mask | (power > self.floor * 100)
        absolute = np.abs(error[evaluated])
        p90 = float(np.quantile(absolute, 0.9)) if absolute.size else 0.0
        maximum = float(absolute.max()) if absolute.size else 0.0
        return dict(
            edges_hz=self.edges.tolist(),
            regions_seconds=self.regions,
            error_db=error.tolist(),
            evaluated=evaluated.tolist(),
            p90_error_db=p90,
            max_error_db=maximum,
            shape_guard=bool(p90 <= 6 and maximum <= 12),
        )
