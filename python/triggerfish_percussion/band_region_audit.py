"""Independent band/time checks, without averaging failures into one score."""

import numpy as np
from scipy.signal import butter, sosfilt

BANDS = (
    (20, 120),
    (120, 250),
    (250, 450),
    (450, 1000),
    (1000, 2000),
    (2000, 4000),
    (4000, 8000),
)
REGIONS = ((0, 0.03), (0.03, 0.1), (0.1, 0.25), (0.25, 0.6), (0.6, 1.2))


class BandRegionAudit:
    """Causal filters with no envelope smoothing across region boundaries.

    Band filters themselves have memory/delay. Reference and synthesis use the same
    filters; these are observable-output measurements, not latent-source estimates.
    """

    def __init__(self, reference, rate):
        self.rate = rate
        self.filters = [
            butter(
                3 if low < 120 else 6,
                (low, min(high, 0.49 * rate)),
                fs=rate,
                btype="bandpass",
                output="sos",
            )
            for low, high in BANDS
        ]
        self.target = self.power(reference)
        self.floor = max(float(self.target.max()), 1e-20) * 1e-8
        # Explicit engineering guard: every audible cell must be within 3 dB.
        # This is necessary, not sufficient, for accepting a timbral fit.
        self.audible = self.target > max(self.floor, 1e-8)

    def power(self, samples):
        result = []
        for sos in self.filters:
            filtered = sosfilt(sos, samples)
            result.append(
                [
                    np.mean(filtered[round(a * self.rate) : round(b * self.rate)] ** 2)
                    for a, b in REGIONS
                ]
            )
        return np.array(result)

    def measure(self, samples):
        power = self.power(samples)
        target_db = 10 * np.log10(np.maximum(self.target, self.floor))
        actual_db = 10 * np.log10(np.maximum(power, self.floor))
        delta = actual_db - target_db
        evaluated = self.audible | (power > max(self.floor, 1e-8))
        return dict(
            bands_hz=BANDS,
            regions_seconds=REGIONS,
            reference_db=target_db.tolist(),
            candidate_db=actual_db.tolist(),
            error_db=delta.tolist(),
            audible=self.audible.tolist(),
            evaluated=evaluated.tolist(),
            worst_audible_error_db=(
                float(np.max(np.abs(delta[evaluated]))) if evaluated.any() else 0.0
            ),
            within_3db=bool(np.all(np.abs(delta[evaluated]) <= 3)),
        )
