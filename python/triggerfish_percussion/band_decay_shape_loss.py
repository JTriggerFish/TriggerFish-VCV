"""Level-independent decay-shape subproblem, followed by absolute-level review."""

import numpy as np

from .transforms import StftConfig, stft


class BandDecayShapeLoss:
    """Fit damping without using bar loudness to conceal a decay error.

    This deliberately removes ONE constant per band in the measurement, never
    from playback. It is not a full-sound acceptance loss. Use an absolute-level
    objective, attack inspection and repeated-hit checks after this subproblem.
    Terminal reference power is a conservative contamination estimate, not a
    claim that the recording's remaining tail is stationary background noise.
    """

    units = "dB shape error"
    influence_threshold = 0.01

    def __init__(self, reference, rate, bands=None):
        if not np.isfinite(rate) or rate <= 0:
            raise ValueError("Expected a positive finite sample rate")
        self.rate, self.frames = rate, len(reference)
        if self.frames / rate < 3.1:
            raise ValueError("Decay-shape fitting needs at least 3.1 seconds")
        self.config = StftConfig(4096, 512)
        default_bands = (
            (100, 300),
            (300, 700),
            (700, 1500),
            (1500, 3000),
            (3000, 6000),
            (6000, 16000),
        )
        self.bands = default_bands if bands is None else tuple(bands)
        if not self.bands or any(
            not np.isfinite((low, high)).all() or not 0 < low < high
            for low, high in self.bands
        ):
            raise ValueError("Expected finite positive ordered frequency bands")
        transformed = stft(reference, rate, self.config)
        self.times = transformed.times_seconds
        self.bin_masks = [
            (transformed.frequencies_hz >= a)
            & (transformed.frequencies_hz < min(b, 0.49 * rate))
            for a, b in self.bands
        ]
        if any(not mask.any() for mask in self.bin_masks):
            raise ValueError("Sample rate cannot resolve the requested bands")
        target = self.power(reference)
        self.floor = max(float(target.max()), 1e-20) * 1e-8
        self.anchor = (self.times >= 0.2) & (self.times < 0.5)
        terminal = (self.times >= self.frames / rate - 0.6) & (
            self.times < self.frames / rate - 0.1
        )
        contamination = np.median(target[:, terminal], axis=1)
        thresholds = np.maximum(target.max(axis=1) * 1e-4, contamination * 10)
        self.mask = (target > thresholds[:, None]) & (
            (self.times >= 0.2) & (self.times < self.frames / rate - 0.1)
        )[None, :]
        if not np.any(self.mask.sum(axis=1) >= 4):
            raise ValueError("Reference has insufficient uncontaminated decay")
        self.target = self.relative_db(target)
        self.specification = dict(
            version="band-decay-shape-v1",
            bands_hz=self.bands,
            fft=4096,
            hop=512,
            anchor_seconds=[0.2, 0.5],
            reference_relative_floor_db=-40,
            terminal_margin_db=10,
            playback_normalization=False,
        )

    def power(self, samples):
        if np.shape(samples) != (self.frames,):
            raise ValueError("Candidate must match the reference duration")
        value = stft(samples, self.rate, self.config)
        return np.array([value.power[mask].sum(axis=0) for mask in self.bin_masks])

    def relative_db(self, power):
        db = 10 * np.log10(np.maximum(power, self.floor))
        return db - db[:, self.anchor].mean(axis=1, keepdims=True)

    def residual(self, samples, regions=range(5)):
        # Search's regional argument is irrelevant to this whole-decay stage.
        error = self.relative_db(self.power(samples)) - self.target
        parts = [
            row[mask] / np.sqrt(mask.sum())
            for row, mask in zip(error, self.mask)
            if mask.sum() >= 4
        ]
        if not parts:
            raise ValueError("Reference has insufficient uncontaminated decay")
        return np.concatenate(parts) / np.sqrt(len(parts))

    def diagnostics(self, samples):
        return dict(shape_error_db=float(np.linalg.norm(self.residual(samples))))
