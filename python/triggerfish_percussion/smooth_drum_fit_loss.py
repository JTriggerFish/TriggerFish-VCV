"""Smooth spectral objective inspired by Schwaer & Mueller (IEEE SPL, 2023).

Flattop windows reduce sidelobe-driven frequency minima; log1p avoids singular
log valleys. Reference-only compression scales preserve level differences.
The existing Hann/dB measurements remain independent audit diagnostics.
"""

import numpy as np
from scipy.signal import ShortTimeFFT, windows

from .short_drum_fit_loss import ShortDrumLoss, REGIONS, REGION_WEIGHTS


class SmoothDrumLoss:
    units = "compressed magnitude"
    influence_threshold = 0.0005
    specification = dict(
        version="smooth-drum-v1",
        window="flattop",
        sizes=[509, 2039, 8191, 16381],
        hop=256,
        compression="log1p(100*magnitude/reference_peak)",
        regions=REGIONS,
        region_weights=REGION_WEIGHTS.tolist(),
    )

    def __init__(self, reference, sample_rate):
        self.audit = ShortDrumLoss(reference, sample_rate)
        self.transforms, self.targets, self.scales, self.times, self.masks = (
            [],
            [],
            [],
            [],
            [],
        )
        for size, maximum in ((509, 16000), (2039, 3000), (8191, 500), (16381, 200)):
            transform = ShortTimeFFT(
                windows.flattop(size, sym=False),
                256,
                sample_rate,
                fft_mode="onesided",
                scale_to="magnitude",
            )
            times = transform.t(len(reference))
            time = (times >= 0) & (times <= len(reference) / sample_rate)
            mask = transform.f <= maximum
            target = np.abs(transform.stft(reference)[mask][:, time])
            scale = 100 / max(float(target.max()), 1e-12)
            self.transforms.append(transform)
            self.targets.append(np.log1p(scale * target))
            self.scales.append(scale)
            self.times.append(times[time])
            self.masks.append((mask, time))

    def residual(self, samples, regions=range(5)):
        result = []
        for transform, target, scale, times, (mask, time) in zip(
            self.transforms, self.targets, self.scales, self.times, self.masks
        ):
            magnitude = np.abs(transform.stft(samples)[mask][:, time])
            error = np.log1p(scale * magnitude) - target
            for region in regions:
                start, end = REGIONS[region]
                selected = (times >= start) & (times < end)
                block = error[:, selected]
                result.append(
                    block.ravel() * np.sqrt(REGION_WEIGHTS[region] / (4 * block.size))
                )
        return np.concatenate(result)

    def diagnostics(self, samples):
        result = self.audit.diagnostics(samples)
        result["smooth_objective"] = float(np.linalg.norm(self.residual(samples)))
        return result
