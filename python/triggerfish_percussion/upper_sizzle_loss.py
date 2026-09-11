"""Regional upper spectrum and modulation, guarded by the accepted low body.

This is a diagnostic composite, not a validated perceptual-equivalence metric.
Reference levels are fixed. Texture normalization never changes audition audio.
"""

import numpy as np
from scipy.signal import stft

from .modal_texture_loss import ModalTextureLoss
from .spectral_bloom_loss import SpectralBloomLoss


class UpperSizzleLoss:
    units = "upper spectral/modulation dB composite with low-body guards"

    def __init__(self, reference, rate):
        if not np.isfinite(rate) or rate < 32000 or len(reference) < 6 * rate:
            raise ValueError("Upper sizzle analysis requires >= 32 kHz and six seconds")
        self.rate, self.frames = rate, len(reference)
        # Equal-ERB-width bins resolve colour changes hidden by 7–14 kHz pooling.
        erb = lambda f: np.log1p(0.00437 * f)
        self.edges = np.expm1(np.linspace(erb(3000), erb(15000), 19)) / 0.00437
        self.regions = [(0.1, 0.35), (0.35, 0.6), (0.6, 1.6), (1.6, 3), (3, 5)]
        self.texture = ModalTextureLoss(
            reference, rate, centres=[3500, 4500, 6000, 8000, 10000, 12000]
        )
        self.bloom = SpectralBloomLoss(reference, rate)
        self.reference = self.measure(reference)
        self.specification = dict(
            version="upper-sizzle-v1",
            edges_hz=self.edges.tolist(),
            regions_seconds=self.regions,
            texture=self.texture.specification,
            selected_modulation_bands_hz=[[8, 32], [32, 128]],
            low_guard_regions_seconds=[(0, 0.05), (0.05, 0.1), *self.regions],
            weights=dict(spectrum=1, texture=0.6, centroid=12, bloom_rise=0.35),
            body_guard_db=2,
            body_guard_weight=8,
            normalization=False,
            search_rejection="Any low-body region changing more than 2 dB on any training seed",
            seed_aggregation="mean plus half the worst seed score",
            interpretation="Search ranking only; audit each metric and each seed",
        )

    def measure(self, audio):
        audio = np.asarray(audio)
        if audio.shape != (self.frames,) or not np.isfinite(audio).all():
            raise ValueError("Expected finite reference-length mono audio")
        f, t, z = stft(
            audio, self.rate, nperseg=4096, noverlap=4096 - round(0.01 * self.rate)
        )
        power = abs(z) ** 2
        selections = [(t >= a) & (t < b) for a, b in self.regions]
        band_power = np.array(
            [
                power[(f >= lo) & (f < hi)].sum(axis=0)
                for lo, hi in zip(self.edges[:-1], self.edges[1:])
            ]
        )
        db = 10 * np.log10(
            np.maximum(
                1e-10, np.array([band_power[:, s].mean(axis=1) for s in selections])
            )
        )
        high = (f >= 3000) & (f <= 15000)
        centroid = []
        for s in selections:
            p = power[high][:, s].mean(axis=1)
            centroid.append(float((f[high] * p).sum() / max(p.sum(), 1e-20)))
        texture, _ = self.texture.features(audio)
        # Middle/late texture, 8–32 and 32–128 Hz; slow bloom is scored separately.
        modulation = 10 * texture.reshape(6, 3, 4)[:, 1:, 1:3]
        low = np.array(
            [
                power[(f >= lo) & (f < hi)].sum(axis=0)
                for lo, hi in [(80, 300), (300, 900)]
            ]
        )
        low_selections = [
            (t >= a) & (t < b) for a, b in [(0, 0.05), (0.05, 0.1), *self.regions]
        ]
        low_db = 10 * np.log10(
            np.maximum(
                1e-10, np.array([low[:, s].mean(axis=1) for s in low_selections])
            )
        )
        return dict(
            spectrum=db.tolist(),
            centroid=centroid,
            modulation=modulation.tolist(),
            low_db=low_db.tolist(),
            rise=self.bloom.diagnostics(audio)["rise_rms_db"],
        )

    def diagnostics(self, audio, baseline=None):
        value = self.measure(audio)
        rms = lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))
        spectral = rms(np.array(value["spectrum"]) - self.reference["spectrum"])
        modulation = rms(np.array(value["modulation"]) - self.reference["modulation"])
        centroid = rms(
            np.log2(
                np.maximum(value["centroid"], 1)
                / np.maximum(self.reference["centroid"], 1)
            )
        )
        low = rms(np.array(value["low_db"]) - baseline["low_db"]) if baseline else 0
        score = (
            spectral
            + 0.6 * modulation
            + 12 * centroid
            + 0.35 * value["rise"]
            + 8 * max(0, low - 2)
        )
        low_peak = (
            float(np.max(np.abs(np.array(value["low_db"]) - baseline["low_db"])))
            if baseline
            else 0
        )
        return dict(
            score=score,
            spectrum_db=spectral,
            modulation_db=modulation,
            centroid_octaves=centroid,
            low_change_db=low,
            rise_db=value["rise"],
            low_max_change_db=low_peak,
            centroids_hz=value["centroid"],
        )

    def score(self, audio):
        return self.diagnostics(audio)["score"]
