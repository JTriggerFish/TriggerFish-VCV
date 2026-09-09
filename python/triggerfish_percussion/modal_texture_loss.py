"""Regional auditory-band envelope modulation, not exact upper ridge positions.

Inspired by auditory texture statistics, not a reproduction of a published
perceptual model. Use alongside fixed-level attack/bloom/decay measurements.
Gain-invariant texture statistics alone cannot judge a percussion fit.
"""

import numpy as np
from scipy.fft import fft, ifft, next_fast_len, fftfreq
from scipy.signal import resample_poly, welch


class ModalTextureLoss:
    units = "regional modulation feature distance"

    def __init__(self, reference, rate, centres=None):
        if not np.isfinite(rate) or rate < 4000 or len(reference) < 2 * rate:
            raise ValueError(
                "Texture analysis requires rate >= 4 kHz and >= two seconds"
            )
        self.rate, self.frames = rate, len(reference)
        self.centres = (
            np.geomspace(1500, min(12000, rate * 0.35), 10)
            if centres is None
            else np.asarray(centres)
        )
        if (
            self.centres.ndim != 1
            or not len(self.centres)
            or not np.isfinite(self.centres).all()
            or np.any(self.centres <= 0)
            or np.any(self.centres >= rate * 0.5)
        ):
            raise ValueError("Texture band centres must be finite and inside Nyquist")
        self.regions = ((0.1, 0.6), (0.6, 1.6), (1.6, min(5.0, len(reference) / rate)))
        self.size = next_fast_len(self.frames * 2)
        self.frequencies = fftfreq(self.size, 1 / rate)
        self.step = max(1, round(rate / 512))
        self.envelope_rate = rate / self.step
        self.masks = []
        for centre in self.centres:
            width = 24.7 * (1 + 0.00437 * centre)
            distance = np.abs(self.frequencies - centre) / width
            self.masks.append(
                np.where(distance < 1, 1 + np.cos(np.pi * np.minimum(distance, 1)), 0)
            )
        self.target, self.power = self.features(reference)
        self.active = self.power > max(self.power.max() * 1e-5, 1e-24)
        if not np.any(self.active):
            raise ValueError(
                "Reference has no measurable texture in the selected bands"
            )
        self.specification = dict(
            version="modal-texture-v1",
            centres_hz=self.centres.tolist(),
            regions=self.regions,
            modulation_bands_hz=[2, 8, 32, 128],
            normalization="band-envelope mean; pair with absolute level loss",
            features=[
                "modulation power 2-8/8-32/32-128 Hz",
                "modulation concentration",
            ],
        )

    def features(self, audio):
        audio = np.asarray(audio)
        if audio.shape != (self.frames,) or not np.isfinite(audio).all():
            raise ValueError(
                "Texture input must be finite mono audio of reference length"
            )
        spectrum = fft(audio, self.size)
        rows, powers = [], []
        for mask in self.masks:
            envelope = np.abs(ifft(spectrum * mask)[: self.frames])
            envelope = resample_poly(envelope, 1, self.step)
            for a, b in self.regions:
                segment = envelope[
                    round(a * self.envelope_rate) : round(b * self.envelope_rate)
                ]
                if len(segment) < 16:
                    raise ValueError(
                        "Texture measurement requires at least two seconds"
                    )
                mean = max(float(np.mean(segment)), 1e-15)
                f, p = welch(
                    segment / mean,
                    self.envelope_rate,
                    nperseg=min(len(segment), 512),
                    detrend="linear",
                )
                power = np.array(
                    [
                        p[(f >= lo) & (f < hi)].sum() * (f[1] - f[0])
                        for lo, hi in ((2, 8), (8, 32), (32, 128))
                    ]
                )
                band = p[(f >= 2) & (f < 128)]
                concentration = band.max() / max(band.sum(), 1e-20)
                rows.append([*np.log10(np.maximum(power, 1e-5)), concentration * 3])
                powers.append(mean * mean)
        return np.asarray(rows), np.asarray(powers)

    def residual(self, audio, regions=None):
        features, _ = self.features(audio)
        error = (features - self.target)[self.active].ravel()
        return error / np.sqrt(max(len(error), 1))

    def score(self, audio):
        return float(np.linalg.norm(self.residual(audio)))

    def diagnostics(self, audio):
        features, power = self.features(audio)
        return dict(
            specification=self.specification,
            target=self.target.tolist(),
            candidate=features.tolist(),
            active=self.active.tolist(),
            candidate_band_power=power.tolist(),
            score=self.score(audio),
        )
