"""Absolute-amplitude attack/body/tail comparison for a single short drum hit.

Three STFT resolutions keep attack noise, body pitch and low ringing visible.
Reference-only floors and weights are immutable during a search. This is an
engineering objective, not a perceptual acceptance test or a pitch tracker.
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .transforms import StftConfig, stft
from .power_envelope import smoothed_power

REGIONS = ((0, 0.03), (0.03, 0.1), (0.1, 0.25), (0.25, 0.6), (0.6, 1.2))
REGION_WEIGHTS = np.array([0.3, 0.3, 0.25, 0.1, 0.05])


def power_db(power, floor):
    return 10 * np.log10(np.maximum(power, floor))


def representations(samples, rate):
    """No zero-padding is used to pretend that frequency resolution improves."""
    result = []
    for size, low, high in ((512, 250, 16000), (2048, 20, 3000), (8192, 0, 500)):
        value = stft(samples, rate, StftConfig(size, 256))
        selected = (value.frequencies_hz >= low) & (value.frequencies_hz <= high)
        power = value.power[selected]
        # Mild bin smoothing reduces dependence on stochastic fine structure.
        result.append((gaussian_filter1d(power, 0.6, axis=0), value.times_seconds))
    # Smooth POWER over 12 ms, sample every 2 ms. A 2 ms power window tracks
    # individual bass cycles and wrongly makes envelope matching phase-sensitive.
    envelope = smoothed_power(samples, 0.012 * rate)[:: round(0.002 * rate)]
    result.append(
        (envelope[None, :], np.arange(len(envelope)) * round(0.002 * rate) / rate)
    )
    return result


class ShortDrumLoss:
    specification = dict(
        version="short-drum-v1",
        windows=[512, 2048, 8192],
        hop=256,
        frequency_ranges=[[250, 16000], [20, 3000], [0, 500]],
        power_smoothing_seconds=0.012,
        power_sample_seconds=0.002,
        floor_db=-70,
        regions=REGIONS,
        region_weights=REGION_WEIGHTS.tolist(),
    )

    def __init__(self, reference, sample_rate):
        self.sample_rate = sample_rate
        self.target = representations(reference, sample_rate)
        self.blocks = []
        for power, times in self.target:
            floor = max(float(power.max()), 1e-20) * 1e-7
            target = power_db(power, floor)
            weight = np.clip((target - target.max() + 50) / 20, 0.1, 1)
            self.blocks.append((target, weight, floor, times))

    def residual(self, samples, regions=range(5)):
        values = representations(samples, self.sample_rate)
        residuals = []
        for (power, _), (target, weight, floor, times) in zip(values, self.blocks):
            error = power_db(power, floor) - target
            for region in regions:
                start, end = REGIONS[region]
                selected = (times >= start) & (times < end)
                weights = weight[:, selected]
                scale = REGION_WEIGHTS[region] / len(values)
                residuals.append(
                    (
                        error[:, selected] * np.sqrt(scale * weights / weights.sum())
                    ).ravel()
                )
        return np.concatenate(residuals)

    def diagnostics(self, samples):
        values = representations(samples, self.sample_rate)
        rows = []
        names = ("attack_spectrum", "body_spectrum", "low_spectrum", "envelope")
        for region, (start, end) in enumerate(REGIONS):
            row = dict(seconds=[start, end])
            for name, (power, _), (target, weight, floor, times) in zip(
                names, values, self.blocks
            ):
                selected = (times >= start) & (times < end)
                error = power_db(power, floor)[:, selected] - target[:, selected]
                row[f"{name}_rmse_db"] = float(
                    np.sqrt(np.average(error**2, weights=weight[:, selected]))
                )
            rows.append(row)
        return dict(
            regions=rows, rms_error_db=float(np.linalg.norm(self.residual(samples)))
        )
