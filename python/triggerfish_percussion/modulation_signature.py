"""Depth, periodicity and cross-band synchrony diagnostics for ringing tails.

These describe modulation, not overall sound quality. Pair with fixed-level
spectral/attack/decay checks; a quieter or noisier signal must not win by itself.
"""

from math import gcd
import numpy as np
from scipy.signal import butter, hilbert, periodogram, resample_poly, sosfiltfilt

BANDS = (
    (90, 180),
    (180, 320),
    (320, 550),
    (550, 900),
    (900, 1600),
    (1600, 2800),
    (2800, 4800),
    (4800, 7000),
)


def modulation_signature(audio, rate):
    """Six seconds mono; inspect 0.5–4.5 s after removing log-linear decay."""
    audio = np.asarray(audio)
    if (
        audio.ndim != 1
        or not np.isfinite(audio).all()
        or not np.isfinite(rate)
        or rate < 16000
    ):
        raise ValueError("Expected finite mono audio at >=16 kHz")
    if rate != int(rate) or len(audio) < 6 * rate:
        raise ValueError("Expected six seconds and an integer sample rate")
    divisor = gcd(int(rate), 16000)
    signal = resample_poly(
        audio[: 6 * int(rate)], 16000 // divisor, int(rate) // divisor
    )
    rows, motion, lines = [], [], []
    for low, high in BANDS:
        band = sosfiltfilt(
            butter(4, [low, high], fs=16000, btype="bandpass", output="sos"), signal
        )
        envelope = np.abs(hilbert(band, 2 * len(band))[: len(band)])
        envelope = resample_poly(envelope, 1, 80)[100:900]
        time = np.arange(len(envelope)) / 200
        trend = np.exp(
            np.polyval(np.polyfit(time, np.log(np.maximum(envelope, 1e-15)), 1), time)
        )
        relative = envelope / np.maximum(trend, 1e-15)
        relative = (
            relative / max(float(relative.mean()), 1e-15) - 1
            if envelope.max() > 1e-15
            else np.zeros_like(envelope)
        )
        f, p = periodogram(relative, fs=200, window="hann", detrend=False)
        select = (f >= 0.5) & (f < 12)
        power = p[select]
        flutter = (f >= 12) & (f < 80)
        flutter_power = p[flutter]
        line_power = p[(f >= 1) & (f < 12)]
        lines.append(line_power / max(float(line_power.sum()), 1e-20))
        motion.append(
            sosfiltfilt(
                butter(3, [0.5, 12], fs=200, btype="bandpass", output="sos"), relative
            )
        )
        rows.append(
            dict(
                band_hz=[low, high],
                depth=float(np.sqrt(power.sum() * (f[1] - f[0]))),
                slow_depth=float(
                    np.sqrt(p[(f >= 0.5) & (f < 3)].sum() * (f[1] - f[0]))
                ),
                fast_depth=float(np.sqrt(p[(f >= 3) & (f < 12)].sum() * (f[1] - f[0]))),
                flutter_depth=float(np.sqrt(flutter_power.sum() * (f[1] - f[0]))),
                flutter_periodicity=float(
                    flutter_power.max() / max(flutter_power.sum(), 1e-20)
                ),
                flutter_dominant_hz=float(f[flutter][np.argmax(flutter_power)]),
                periodicity=float(power.max() / max(power.sum(), 1e-20)),
                dominant_hz=float(f[select][np.argmax(power)]),
                level=float(np.sqrt(np.mean(band[8000:72000] ** 2))),
            )
        )
    active = [
        i
        for i, r in enumerate(rows)
        if r["level"] > max(x["level"] for x in rows) * 0.003
        and np.std(motion[i]) > 1e-9
    ]
    correlations = [
        float(abs(np.corrcoef(motion[i], motion[j])[0, 1]))
        for i in active
        for j in active
        if i < j
    ]
    shared = float(np.mean([lines[i] for i in active], axis=0).max()) if active else 0
    return dict(
        bands=rows,
        common_line_strength=shared,
        mean_synchrony=float(np.mean(correlations)) if correlations else 0,
        max_synchrony=max(correlations, default=0),
        ranges_hz=dict(slow=[0.5, 3], fast=[3, 12], flutter=[12, 80]),
        interpretation="Descriptive only; common bloom can also correlate bands",
    )


def excess_motion(candidate, reference):
    """Penalize excess depth/regularity, never reward making the tail silent."""
    active = [
        i
        for i, r in enumerate(reference["bands"])
        if r["level"] > max(b["level"] for b in reference["bands"]) * 0.003
    ]
    if not active:
        raise ValueError("Reference has no measurable modulation bands")
    depth = [
        max(
            0,
            np.log2(
                max(candidate["bands"][i]["depth"], 0.001)
                / max(reference["bands"][i]["depth"], 0.001)
            ),
        )
        for i in active
    ]
    periodicity = [
        max(
            0,
            candidate["bands"][i]["periodicity"] - reference["bands"][i]["periodicity"],
        )
        for i in active
    ]
    return float(
        np.mean(depth)
        + 2 * np.mean(periodicity)
        + 2
        * max(0, candidate["common_line_strength"] - reference["common_line_strength"])
    )
