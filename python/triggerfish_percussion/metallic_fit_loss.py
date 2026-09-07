"""Long metallic trajectories plus explicit attack and upper-ringing constraints."""

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .trajectory_fit_loss import REGIONS, TrajectoryLoss, to_db
from .transforms import StftConfig, stft


class MetallicFitLoss:
    """Fixed-reference dB residuals, never a listening-acceptance criterion.

    The long trajectory alone blurs the contact and ignores narrow upper ridges.
    Add 512-sample attack frames and 4096-sample regional spectra through 16 kHz.
    Slight spectral smoothing tolerates seed beating without erasing ringing.
    """

    def __init__(self, reference, sample_rate):
        self.rate = sample_rate
        self.trajectory = TrajectoryLoss(reference, sample_rate)
        self.target = self.extra_features(reference)
        self.floors = [max(float(x.max()), 1e-20) * 1e-7 for x in self.target]
        self.db = [to_db(x, f) for x, f in zip(self.target, self.floors)]
        self.weights = [np.clip((x - x.max() + 55) / 25, 0.02, 1) for x in self.db]
        self.specification = dict(
            version="metallic-attack-trajectory-ridges-v1",
            regions=REGIONS,
            shares=dict(trajectory=0.6, attack=0.2, upper_ridges=0.2),
            attack_window=512,
            attack_hop=128,
            attack_seconds=0.12,
            ridge_window=4096,
            ridge_hop=1024,
            ridge_smoothing_bins=0.8,
            reference_floor_db=-70,
            normalization=False,
        )

    def extra_features(self, samples):
        attack = stft(
            samples[: round(0.14 * self.rate)], self.rate, StftConfig(512, 128)
        )
        selected = (attack.frequencies_hz >= 80) & (attack.frequencies_hz <= 16000)
        onset = attack.power[selected][:, attack.times_seconds < 0.12]
        long = stft(samples, self.rate, StftConfig(4096, 1024))
        selected = (long.frequencies_hz >= 1000) & (long.frequencies_hz <= 16000)
        ridges = np.array(
            [
                long.power[selected][
                    :, (long.times_seconds >= a) & (long.times_seconds < b)
                ].mean(axis=1)
                for a, b in REGIONS
            ]
        )
        return onset, gaussian_filter1d(ridges, 0.8, axis=1)

    def residual(self, samples, regions=range(5)):
        regions = tuple(regions)
        parts = [np.sqrt(0.6) * self.trajectory.residual(samples, regions)]
        features = self.extra_features(samples)
        if 0 in regions:
            w = self.weights[0]
            parts.append(
                (
                    np.sqrt(0.2 * w / w.sum())
                    * (to_db(features[0], self.floors[0]) - self.db[0])
                ).ravel()
            )
        w = self.weights[1][list(regions)]
        error = (
            to_db(features[1], self.floors[1])[list(regions)]
            - self.db[1][list(regions)]
        )
        parts.append((np.sqrt(0.2 * w / w.sum()) * error).ravel())
        return np.concatenate(parts)

    def diagnostics(self, samples):
        result = self.trajectory.diagnostics(samples)
        result["trajectory_rms_db"] = result.pop("rms_error_db")
        result["rms_error_db"] = float(np.linalg.norm(self.residual(samples)))
        return result
