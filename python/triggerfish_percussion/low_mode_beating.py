"""Slow low-band beating after removal of a fitted exponential decay.

This diagnostic distinguishes slow pulses from fast flutter; it is not a
complete perceptual loss. Compare absolute spectra/decays separately. Four
seconds gives 0.25-Hz modulation resolution, unlike the short upper-wash test.
"""

import numpy as np
from scipy.signal import butter, hilbert, periodogram, resample_poly, sosfiltfilt
from math import gcd


class LowModeBeating:
    bands = ((90, 180), (180, 320), (320, 550), (550, 900))
    modulation_bands = ((0.5, 3), (3, 8), (8, 30))

    def __init__(self, reference, rate, region=(0.5, 4.5)):
        self.rate, self.frames, self.region = rate, len(reference), region
        if (
            not np.isfinite(rate)
            or rate != int(rate)
            or rate < 4000
            or not 0 <= region[0] < region[1] <= len(reference) / rate
        ):
            raise ValueError("Expected integer sample rate >= 4 kHz and valid region")
        if region[1] - region[0] < 4:
            raise ValueError("Slow beating requires at least four seconds")
        self.target = self.analyze(reference)
        if max(r["band_rms"] for r in self.target) < 1e-12:
            raise ValueError("Reference has no measurable low-band energy")
        self.specification = dict(
            version="low-mode-beating-v1",
            bands_hz=self.bands,
            region_seconds=region,
            modulation_bands_hz=self.modulation_bands,
            detrend="log-envelope linear regression (exponential decay)",
            modulation_resolution_hz=1 / (region[1] - region[0]),
        )

    def analyze(self, audio):
        audio = np.asarray(audio)
        if audio.shape != (self.frames,) or not np.isfinite(audio).all():
            raise ValueError("Expected finite mono audio of reference length")
        divisor = gcd(int(self.rate), 4000)
        audio = resample_poly(audio, 4000 // divisor, int(self.rate) // divisor)
        rows = []
        for low, high in self.bands:
            filtered = sosfiltfilt(
                butter(4, [low, high], fs=4000, btype="bandpass", output="sos"), audio
            )
            # Padding avoids joining the decayed tail to the attack in Hilbert.
            envelope = np.abs(hilbert(filtered, 2 * len(filtered))[: len(filtered)])
            envelope = resample_poly(envelope, 1, 20)
            a, b = [round(t * 200) for t in self.region]
            segment = np.maximum(envelope[a:b], 1e-15)
            time = np.arange(len(segment)) / 200
            trend = np.exp(np.polyval(np.polyfit(time, np.log(segment), 1), time))
            relative = segment / trend
            # Normalize modulation depth, not reference/candidate audio gain.
            relative = relative / relative.mean() - 1
            f, p = periodogram(relative, fs=200, window="hann", detrend=False)
            df = f[1] - f[0]
            power = [
                float(p[(f >= lo) & (f < hi)].sum() * df)
                for lo, hi in self.modulation_bands
            ]
            selected = (f >= 0.5) & (f < 30)
            rows.append(
                dict(
                    band_hz=[low, high],
                    power=power,
                    dominant_hz=float(f[selected][np.argmax(p[selected])]),
                    fast_fraction=sum(power[1:]) / max(sum(power), 1e-20),
                    relative_envelope=relative.tolist(),
                    spectrum=p.tolist(),
                    band_rms=float(
                        np.sqrt(
                            np.mean(
                                filtered[
                                    round(self.region[0] * 4000) : round(
                                        self.region[1] * 4000
                                    )
                                ]
                                ** 2
                            )
                        )
                    ),
                )
            )
        return rows

    def score_rows(self, rows):
        target = np.array([r["power"] for r in self.target])
        actual = np.array([r["power"] for r in rows])
        active = np.array([r["band_rms"] for r in self.target])
        active = active > max(active.max() * 1e-3, 1e-12)
        if not active.any():
            raise ValueError("Reference has no measurable low-band energy")
        error = np.log10(np.maximum(actual[active], 1e-5)) - np.log10(
            np.maximum(target[active], 1e-5)
        )
        return float(np.sqrt(np.mean(error**2)))

    def score(self, audio):
        return self.score_rows(self.analyze(audio))
