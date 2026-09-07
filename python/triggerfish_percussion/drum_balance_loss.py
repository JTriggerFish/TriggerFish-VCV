"""Broad-band envelope constraints complement, not replace, spectral matching.

Identical causal filters measure reference and candidate. Their group delay is
part of the measurement, not interpreted as an instrument's excitation delay.
One reference-derived floor keeps inaudible upper-band tails from dominating.
"""

import numpy as np
from scipy.signal import butter, sosfilt

from .power_envelope import smoothed_power
from .short_drum_fit_loss import ShortDrumLoss, REGIONS, REGION_WEIGHTS, power_db

BANDS = (
    (20, 45),
    (45, 90),
    (90, 180),
    (180, 350),
    (350, 700),
    (700, 1400),
    (1400, 2800),
    (2800, 5600),
    (5600, 16000),
)


class DrumBalanceLoss:
    """Equal shares of the existing spectral objective and band envelopes."""

    def __init__(self, reference, sample_rate):
        self.spectral = ShortDrumLoss(reference, sample_rate)
        self.rate = sample_rate
        self.hop = max(1, round(0.004 * sample_rate))
        self.bands = [
            (lo, min(hi, 0.49 * sample_rate))
            for lo, hi in BANDS
            if lo < 0.49 * sample_rate
        ]
        self.filters = [
            butter(2, band, btype="bandpass", fs=sample_rate, output="sos")
            for band in self.bands
        ]
        self.target_power = self.envelopes(reference)
        self.floor = max(float(self.target_power.max()), 1e-20) * 1e-6
        self.target = power_db(self.target_power, self.floor)
        self.weight = np.clip((self.target - self.target.max() + 50) / 30, 0.05, 1)
        self.times = np.arange(self.target.shape[1]) * self.hop / sample_rate
        self.specification = dict(
            version="drum-balance-v1",
            spectral=self.spectral.specification,
            bands_hz=self.bands,
            butterworth_order=2,
            causal=True,
            smoothing_seconds=dict(bass=0.02, other=0.012),
            hop_seconds=self.hop / sample_rate,
            common_reference_floor_db=-60,
            spectral_weight=0.5,
            envelope_weight=0.5,
        )

    def envelopes(self, samples):
        return np.array(
            [
                smoothed_power(
                    sosfilt(sos, samples), (0.02 if band[0] < 90 else 0.012) * self.rate
                )[:: self.hop]
                for band, sos in zip(self.bands, self.filters)
            ]
        )

    def balance_residual(self, samples, regions=range(5)):
        error = power_db(self.envelopes(samples), self.floor) - self.target
        residuals = []
        for region in regions:
            start, end = REGIONS[region]
            selected = (self.times >= start) & (self.times < end)
            weight = self.weight[:, selected]
            residuals.append(
                (
                    error[:, selected]
                    * np.sqrt(REGION_WEIGHTS[region] * weight / weight.sum())
                ).ravel()
            )
        return np.concatenate(residuals)

    def residual(self, samples, regions=range(5)):
        regions = tuple(regions)
        return np.concatenate(
            (
                self.spectral.residual(samples, regions),
                self.balance_residual(samples, regions),
            )
        ) / np.sqrt(2)

    def diagnostics(self, samples):
        result = self.spectral.diagnostics(samples)
        result["spectral_rms_error_db"] = result.pop("rms_error_db")
        result["rms_error_db"] = float(np.linalg.norm(self.residual(samples)))
        result["band_envelope_rmse_db"] = float(
            np.linalg.norm(self.balance_residual(samples))
        )
        power = self.envelopes(samples)
        rows = []
        for region in REGIONS:
            selected = (self.times >= region[0]) & (self.times < region[1])
            reference = power_db(
                self.target_power[:, selected].mean(axis=1), self.floor
            )
            candidate = power_db(power[:, selected].mean(axis=1), self.floor)
            rows.append(
                dict(
                    seconds=region,
                    reference_db=reference.tolist(),
                    candidate_db=candidate.tolist(),
                    difference_db=(candidate - reference).tolist(),
                )
            )
        result["band_envelopes"] = dict(bands_hz=self.bands, regions=rows)
        return result
