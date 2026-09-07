"""Search residuals paired with independent, non-averaged output checks."""

import numpy as np
from .region_spectrum_audit import RegionSpectrumAudit


class RegionFitLoss:
    units = "band/3dB and shape/6dB"
    influence_threshold = 0  # retain weak new handles; measure every direction

    def __init__(self, audit, reference, rate):
        self.audit = audit
        self.shape = RegionSpectrumAudit(reference, rate)
        self.specification = dict(
            version="region-fit-v2",
            band_tolerance_db=3,
            shape_tolerance_db=6,
            normalization="none",
            quiet_regions="one-sided excess; 3dB margin below evaluation threshold",
        )

    @staticmethod
    def errors(power, target, floor, active, threshold):
        power = np.maximum(power, floor)
        signed = 10 * np.log10(power / np.maximum(target, floor))
        # Fixed-size residuals also penalize new energy outside the reference's
        # active cells. No candidate-dependent masking or silence exclusion.
        ceiling = np.maximum(target, max(floor, threshold * 0.5))
        excess = np.maximum(0, 10 * np.log10(power / ceiling))
        return np.where(active, signed, excess)

    def residual(self, samples, regions=range(5)):
        if tuple(regions) != tuple(range(5)):
            raise ValueError("RegionFitLoss requires all regions")
        band = self.errors(
            self.audit.power(samples),
            self.audit.target,
            self.audit.floor,
            self.audit.audible,
            max(self.audit.floor, 1e-8),
        )
        shape = self.errors(
            self.shape.power(samples),
            self.shape.target,
            self.shape.floor,
            self.shape.mask,
            self.shape.floor * 100,
        )
        return np.concatenate((band.ravel() / 3, shape.ravel() / 6))

    def diagnostics(self, samples):
        return dict(band=self.audit.measure(samples), shape=self.shape.measure(samples))
