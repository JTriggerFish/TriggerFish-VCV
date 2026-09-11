"""Explicit sound-design target for the gong texture listening experiments.

Retain the incoming low body; compare upper texture to the real recording.
This is not an assertion that the incoming bass matches that recording.
"""

import numpy as np

from triggerfish_percussion.regional_spectrum_audit import RegionalSpectrumAudit
from audit_gong_sizzle import ridge_contrast


class GongTextureComparison:
    units = "regional spectral RMS dB; explicit hybrid sound-design target"
    regions = ((0.04, 0.2), (0.2, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 6))

    def __init__(self, reference, body, rate):
        self.rate = rate
        self.spectrum = RegionalSpectrumAudit(reference, rate, self.regions)
        self.edges = self.spectrum.edges
        self.centres = np.sqrt(self.edges[:-1] * self.edges[1:])
        mix = np.clip(np.log2(self.centres / 800) / np.log2(1800 / 800), 0, 1)
        mix = mix * mix * (3 - 2 * mix)
        self.target_power = np.exp(
            (1 - mix) * np.log(np.maximum(self.spectrum.power(body), 1e-24))
            + mix * np.log(np.maximum(self.spectrum.power(reference), 1e-24))
        )
        self.floor = (
            np.maximum(self.target_power.max(axis=1, keepdims=True), 1e-20) * 1e-6
        )
        self.target = self.db(self.target_power)
        self.weights = np.maximum(
            0.05,
            np.sqrt(self.target_power / self.target_power.max(axis=1, keepdims=True)),
        )
        self.target_ridges = np.array(ridge_contrast(reference, rate))
        self.specification = dict(
            purpose="sound-design refinement, not full reference calibration",
            body="incoming gong with output EQ bypassed",
            upper="unaltered recording at saved reference gain",
            transition_hz=[800, 1800],
            regions_seconds=self.regions,
            frequency_edges_hz=self.edges.tolist(),
            normalization=False,
            observation="specified by fitting history; spectral diagnostic does not choose search coordinates",
        )

    def db(self, power):
        return 10 * np.log10(np.maximum(power, self.floor))

    def power_residual(self, power):
        return (self.db(power) - self.target) * np.sqrt(self.weights)

    def residual(self, audio):
        return self.power_residual(self.spectrum.power(audio)).ravel()

    def score(self, audio):
        return float(np.sqrt(np.mean(self.residual(audio) ** 2)))

    def diagnostics(self, audio):
        return self.measure(audio)

    def measure(self, audio):
        error = self.db(self.spectrum.power(audio)) - self.target
        high = self.centres >= 1800
        low = self.centres < 800
        weights = self.weights
        bias = (error * weights).sum(axis=0) / weights.sum(axis=0)
        shape = error - bias
        rms = lambda x, w: float(np.sqrt(np.sum(x * x * w) / np.sum(w)))
        ridges = np.array(ridge_contrast(audio, self.rate))
        return dict(
            body_db=rms(error[:, low], weights[:, low]),
            upper_db=rms(error[:, high], weights[:, high]),
            upper_shape_db=rms(shape[:, high], weights[:, high]),
            upper_bias_db=float(np.mean(bias[high])),
            ridge_contrast_db=ridges.tolist(),
            ridge_error_db=float(np.sqrt(np.mean((ridges - self.target_ridges) ** 2))),
        )
