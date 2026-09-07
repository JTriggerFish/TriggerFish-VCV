"""Resolve narrow ringing missed by broad-band kick envelope averages.

Reference and candidate use identical centered STFTs. A fixed reference-derived
floor and frozen audibility mask prevent a candidate from hiding its own peaks.
This is a diagnostic spectral objective, not a perceptual-quality certificate.
"""

import numpy as np
from scipy.ndimage import uniform_filter1d

from .drum_balance_loss import DrumBalanceLoss
from .transforms import StftConfig, stft


class RidgeBalanceLoss:
    influence_threshold = 0.005
    regions = ((0, 0.04), (0.04, 0.08), (0.08, 0.16), (0.16, 0.26), (0.26, 0.4))

    def __init__(self, reference, baseline, rate, low_resolution=False):
        self.rate = rate
        self.balance = DrumBalanceLoss(reference, rate)
        self.config = StftConfig(8192 if low_resolution else 4096, 256)
        value = stft(reference, rate, self.config)
        self.bins = (value.frequencies_hz >= (20 if low_resolution else 100)) & (
            value.frequencies_hz <= (250 if low_resolution else 4000)
        )
        self.times = value.times_seconds
        self.target = self.power(reference)
        self.floor = max(self.target.max(), 1e-20) * 1e-6
        self.target_db = self.db(self.target)
        # Include erroneous baseline peaks even where the reference is quiet.
        self.mask = np.maximum(self.target, self.power(baseline)) > self.floor * 10
        self.weights = np.where(self.mask, 1.0, 0.1)
        self.low = (
            None
            if low_resolution
            else RidgeBalanceLoss(reference, baseline, rate, True)
        )
        self.specification = dict(
            version="ridge-balance-v1",
            fft=4096,
            hop=256,
            regions=self.regions,
            hz=[100, 4000],
            frequency_smoothing_bins=3,
            reference_floor_db=-60,
            weighted_excess_knee=3,
            excess_residual_scale=1,
            bass_fft=8192,
            bass_hz=[20, 250],
            inactive_weight=0.1,
            underlying=self.balance.specification,
        )

    def power(self, samples):
        value = stft(samples, self.rate, self.config)
        return uniform_filter1d(
            value.power[self.bins],
            1 if self.config.window_samples == 8192 else 3,
            axis=0,
            mode="nearest",
        )

    def db(self, power):
        return 10 * np.log10(np.maximum(power, self.floor))

    def ridge_residual(self, samples, regions):
        error = self.db(self.power(samples)) - self.target_db
        result = []
        for region in regions:
            start, end = self.regions[region]
            selected = (self.times >= start) & (self.times < end)
            weights = self.weights[:, selected]
            e = (error[:, selected] * np.sqrt(weights)).ravel()
            # Knee is in weighted residual units: 3 dB in active bins,
            # 3/sqrt(.1) dB in quiet bins. Preserve the fitted v1 objective.
            result.extend(
                (
                    e / np.sqrt(max(1, len(e))),
                    np.maximum(e - 3, 0) / np.sqrt(max(1, len(e))),
                )
            )
        return np.concatenate(result) / np.sqrt(len(tuple(regions)))

    def residual(self, samples, regions=range(5)):
        regions = tuple(regions)
        return np.concatenate(
            (
                self.balance.residual(samples, regions),
                self.ridge_residual(samples, regions),
                self.low.ridge_residual(samples, regions),
            )
        ) / np.sqrt(3)

    def diagnostics(self, samples):
        result = self.balance.diagnostics(samples)
        result["ridge_error"] = float(
            np.linalg.norm(self.ridge_residual(samples, range(5)))
        )
        result["bass_ridge_error"] = float(
            np.linalg.norm(self.low.ridge_residual(samples, range(5)))
        )
        result["combined_error"] = float(np.linalg.norm(self.residual(samples)))
        return result
