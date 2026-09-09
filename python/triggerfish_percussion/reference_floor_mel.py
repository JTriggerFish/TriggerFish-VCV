"""Auraloss Mel comparison with reference-fixed per-resolution power floors.

This bounds the log-magnitude loss's effective dynamic range. It is not an
absolute hearing threshold and neither audio nor display levels are changed.
"""

import numpy as np

from .perceptual_fit_losses import AuralossMel


class ReferenceFloorMel(AuralossMel):
    def __init__(self, reference, rate, dynamic_range_db=60):
        reference = np.asarray(reference)
        if (
            reference.ndim != 1
            or len(reference) <= 4096
            or not np.isfinite(reference).all()
            or not np.any(reference)
            or not np.isfinite(dynamic_range_db)
            or not 30 <= dynamic_range_db <= 100
        ):
            raise ValueError(
                "Expected finite nonzero mono reference and 30–100 dB range"
            )
        super().__init__(reference, rate)
        # Match the offline autograd adapter's precision: FFT cancellation near
        # the floor otherwise introduces a float32/64 objective discrepancy.
        self.loss = self.loss.double()
        self.target = self.torch.tensor(reference, dtype=self.torch.float64)[None, None]
        floors = []
        with self.torch.inference_mode():
            for transform in self.loss.stft_losses:
                transform.eps = (
                    0.0  # Measure actual power, not the library's default floor.
                )
                magnitude, _ = transform.stft(self.target.reshape(1, -1))
                # Library eps acts on linear FFT POWER, before sqrt and Mel.
                # Each window's FFT scale differs; derive its floor separately.
                floor = max(
                    float(magnitude.max()) ** 2 * 10 ** (-dynamic_range_db / 10), 1e-30
                )
                transform.eps = floor
                floors.append(floor)
        self.specification.update(
            comparison="reference-floor-mel-v2",
            scoring_dtype="float64",
            dynamic_range_db=dynamic_range_db,
            floor_domain="linear STFT power before Mel projection",
            eps=floors,
            normalization=False,
        )

    def score(self, samples):
        samples = np.asarray(samples)
        if samples.ndim != 1 or not np.isfinite(samples).all():
            raise ValueError("Expected finite mono candidate")
        if len(samples) != self.target.shape[-1]:
            raise ValueError("Candidate/reference lengths must match")
        with self.torch.inference_mode():
            value = self.torch.tensor(samples, dtype=self.torch.float64)[None, None]
            return float(self.loss(value, self.target))
