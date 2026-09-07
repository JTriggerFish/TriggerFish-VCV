"""Shared metallic objective: energy, loud components, contrast and attack.

Linear magnitude complements log energy (DDSP/auraloss); local log contrast
explicitly tests resolved ringing versus wash. The weights are engineering
choices tested across instruments, not a validated perceptual distance.
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .metallic_balance_features import ATTACK_BINS, MetallicBalanceFeatures
from .trajectory_fit_loss import to_db


class MetallicBalanceLoss:
    units = "composite (dB and scaled relative magnitude)"

    def __init__(
        self, reference, sample_rate, contrast_weighting="linear", fast_attack=False
    ):
        if contrast_weighting not in ("linear", "erb"):
            raise ValueError("Contrast weighting must be linear or erb")
        self.features = MetallicBalanceFeatures(sample_rate, len(reference))
        self.target = self.features(reference)
        self.floors = {
            key: max(float(value.max()), 1e-20) * 1e-7
            for key, value in self.target.items()
        }
        self.db = {
            key: to_db(value, self.floors[key]) for key, value in self.target.items()
        }
        self.weights = {
            key: np.clip((value - value.max() + 60) / 30, 0.02, 1)
            for key, value in self.db.items()
        }
        if contrast_weighting == "erb":
            # Equal ERB-rate intervals, not equal counts of linear FFT bins.
            # d ERB / d f is proportional to 1/(f+1/.00437); the common
            # constant cancels in weighted mean squares.
            self.weights["spectrum"] /= self.features.frequencies[None, :] + 1 / 0.00437
        self.contrast = self.local_contrast(self.db["spectrum"])
        self.specification = dict(
            version="metallic-balance-v2",
            contrast_weighting=contrast_weighting,
            fast_attack=bool(fast_attack),
            regions=self.features.regions,
            shares=dict(envelope=0.35, linear_spectrum=0.3, contrast=0.2, attack=0.15),
            linear_scale=20,
            resolved_window=16384,
            resolved_hop=4096,
            contrast_smoothing_bins=8,
            power_smoothing_seconds=0.012,
            causal_band_filters=True,
            power_smoothing="centered Gaussian",
            bands_hz=self.features.bands,
            reference_floor_db=-70,
            normalization=False,
        )

    @staticmethod
    def local_contrast(db):
        return db - gaussian_filter1d(db, 8, axis=1)

    @staticmethod
    def weighted(error, weight):
        return (error * np.sqrt(weight / weight.sum())).ravel()

    def components(self, samples, regions=range(5)):
        regions = tuple(regions)
        values = self.features(samples)
        db = {key: to_db(value, self.floors[key]) for key, value in values.items()}
        contrast = self.local_contrast(db["spectrum"])
        parts = dict(envelope=[], linear_spectrum=[], contrast=[], attack=[])
        for r in regions:
            a, b = self.features.regions[r]
            selected = (self.features.times >= a) & (self.features.times < b)
            # Split the contact from the rest of the first region: don't let
            # 90 ms of decay swamp a missing 30 ms attack.
            masks = [selected]
            if r == 0:
                masks = [
                    selected & (self.features.times < 0.03),
                    selected & (self.features.times >= 0.03),
                ]
            for mask in masks:
                parts["envelope"].append(
                    self.weighted(
                        db["envelope"][:, mask] - self.db["envelope"][:, mask],
                        self.weights["envelope"][:, mask],
                    )
                    / np.sqrt(len(masks) * len(regions))
                )
            ref = np.sqrt(self.target["spectrum"][r])
            difference = np.sqrt(values["spectrum"][r]) - ref
            # A reference-only denominator scales residual units; it does NOT
            # normalize either waveform or remove candidate gain errors.
            denominator = max(
                np.linalg.norm(ref), np.sqrt(self.floors["spectrum"] * len(ref))
            )
            parts["linear_spectrum"].append(
                20 * difference / denominator / np.sqrt(len(regions))
            )
            parts["contrast"].append(
                self.weighted(
                    contrast[r] - self.contrast[r], self.weights["spectrum"][r]
                )
                / np.sqrt(len(regions))
            )
        if 0 in regions:
            weight = 0.5 if self.specification["fast_attack"] else 1
            parts["attack"].append(
                np.sqrt(weight)
                * self.weighted(
                    db["attack"] - self.db["attack"], self.weights["attack"]
                )
            )
            if self.specification["fast_attack"]:
                parts["attack"].append(
                    np.sqrt(0.5 / len(db["transient"]))
                    * (db["transient"] - self.db["transient"])
                )
        return {key: np.concatenate(value) for key, value in parts.items() if value}

    def residual(self, samples, regions=range(5)):
        parts = self.components(samples, regions)
        return np.concatenate(
            [
                np.sqrt(self.specification["shares"][key]) * value
                for key, value in parts.items()
            ]
        )

    def diagnostics(self, samples):
        samples = np.asarray(samples, dtype=np.float64)
        parts = self.components(samples)
        power = np.array(
            [
                np.mean(
                    samples[
                        round(a * self.features.rate) : round(b * self.features.rate)
                    ]
                    ** 2
                )
                for a, b in ATTACK_BINS
            ]
        )
        return dict(
            rms_error_db=float(
                np.sqrt(
                    sum(
                        self.specification["shares"][key] * np.dot(value, value)
                        for key, value in parts.items()
                    )
                )
            ),
            units=self.units,
            attack_bins_seconds=ATTACK_BINS,
            attack_bin_reference_db=self.db["transient"].tolist(),
            attack_bin_difference_db=(
                to_db(power, self.floors["transient"]) - self.db["transient"]
            ).tolist(),
            components={
                key: float(np.linalg.norm(value)) for key, value in parts.items()
            },
        )
