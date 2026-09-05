"""Reference-anchored, full-duration diagnostics for stochastic percussion.

No candidate normalization, time warp, per-band level matching or phase loss.
Each residual block has explicit units (dB) and equal mean-square weight.
"""

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .transforms import StftConfig, stft

REGIONS = ((0, 0.12), (0.12, 0.5), (0.5, 1.5), (1.5, 3), (3, 6))


def to_db(power, floor):
    return 10 * np.log10(np.maximum(power, floor))


@dataclass
class TrajectoryFeatures:
    bands: np.ndarray
    ridges: np.ndarray
    times: np.ndarray
    frequencies: np.ndarray
    band_centres: np.ndarray


def features(samples, sample_rate):
    short = stft(samples, sample_rate, StftConfig(2048, 512))
    short_power = short.power
    # Uniform ERB-rate edges, with frequencies clipped to actual Nyquist.
    erb = lambda hz: 21.4 * np.log10(1 + 0.00437 * hz)
    edges = (
        10 ** (np.linspace(erb(40), erb(min(16000, sample_rate / 2)), 37) / 21.4) - 1
    ) / 0.00437
    bands = np.array(
        [
            short_power[
                (short.frequencies_hz >= low) & (short.frequencies_hz < high)
            ].sum(axis=0)
            for low, high in zip(edges[:-1], edges[1:])
        ]
    )
    # 32 ms power smoothing suppresses seed-specific beating, not bloom shape.
    bands = gaussian_filter1d(bands, 0.032 * sample_rate / 512, axis=1)
    long = stft(samples, sample_rate, StftConfig(8192, 1024))
    selected = (long.frequencies_hz >= 40) & (long.frequencies_hz <= 3000)
    long_power = long.power[selected]
    ridges = np.array(
        [
            np.mean(
                long_power[
                    :, (long.times_seconds >= start) & (long.times_seconds < end)
                ],
                axis=1,
            )
            for start, end in REGIONS
        ]
    )
    return TrajectoryFeatures(
        bands,
        gaussian_filter1d(ridges, 0.6, axis=1),
        short.times_seconds,
        long.frequencies_hz[selected],
        np.sqrt(edges[:-1] * edges[1:]),
    )


class TrajectoryLoss:
    def __init__(self, reference, sample_rate):
        self.sample_rate = sample_rate
        self.target = features(reference, sample_rate)
        self.floor = max(np.max(self.target.bands), np.max(self.target.ridges)) * 1e-8
        self.band_db = to_db(self.target.bands, self.floor)
        self.ridge_db = to_db(self.target.ridges, self.floor)
        # Reference-only salience weighting; missing candidate energy is penalized.
        self.band_weight = np.clip(
            (self.band_db - self.band_db.max() + 60) / 20, 0.05, 1
        )
        self.ridge_weight = np.clip(
            (self.ridge_db - self.ridge_db.max() + 50) / 20, 0.05, 1
        )

    def residual(self, samples, regions=range(5)):
        value = features(samples, self.sample_rate)
        band_error = to_db(value.bands, self.floor) - self.band_db
        ridge_error = to_db(value.ridges, self.floor) - self.ridge_db
        result = []
        for region in regions:
            start, end = REGIONS[region]
            selected = (self.target.times >= start) & (self.target.times < end)
            weight = self.band_weight[:, selected]
            result.append(
                (band_error[:, selected] * np.sqrt(weight / weight.sum())).ravel()
            )
            weight = self.ridge_weight[region]
            result.append(ridge_error[region] * np.sqrt(weight / weight.sum()))
        return np.concatenate(result) / np.sqrt(len(result))

    def diagnostics(self, samples):
        value = features(samples, self.sample_rate)
        error = to_db(value.bands, self.floor) - self.band_db
        rows = []
        for index, (start, end) in enumerate(REGIONS):
            selected = (self.target.times >= start) & (self.target.times < end)
            weight = self.band_weight[:, selected]
            ridge_error = to_db(value.ridges[index], self.floor) - self.ridge_db[index]
            rows.append(
                dict(
                    seconds=[start, end],
                    band_rmse_db=float(
                        np.sqrt(np.sum(error[:, selected] ** 2 * weight) / weight.sum())
                    ),
                    ridge_rmse_db=float(
                        np.sqrt(
                            np.average(ridge_error**2, weights=self.ridge_weight[index])
                        )
                    ),
                )
            )
        groups = ((40, 300), (300, 1000), (1000, 3000), (3000, 8000), (8000, 16000))
        peaks = lambda item: [
            float(
                item.times[
                    np.argmax(
                        item.bands[
                            (item.band_centres >= low) & (item.band_centres < high)
                        ].sum(axis=0)
                    )
                ]
            )
            for low, high in groups
        ]
        rms = np.sqrt(
            np.mean(
                [
                    row[key] ** 2
                    for row in rows
                    for key in ("band_rmse_db", "ridge_rmse_db")
                ]
            )
        )
        return dict(
            regions=rows,
            rms_error_db=float(rms),
            reference_peak_times=peaks(self.target),
            candidate_peak_times=peaks(value),
        )
