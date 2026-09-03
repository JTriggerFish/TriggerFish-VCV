"""Print compact, level-aligned time/frequency diagnostics for two WAV files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from triggerfish_percussion.alignment import detect_impact_onset
from triggerfish_percussion.audio_io import read_wav, resample_audio

BANDS = ((20, 150), (150, 400), (400, 1000), (1000, 3000), (3000, 8000), (8000, 20000))
REGIONS = ((0.0, 0.015), (0.015, 0.12), (0.12, 0.6), (0.6, 1.5))


def _region(
    signal: np.ndarray, onset: int, rate: int, limits: tuple[float, float]
) -> np.ndarray:
    first = onset + round(limits[0] * rate)
    last = min(signal.size, onset + round(limits[1] * rate))
    return signal[first:last]


def _spectrum(signal: np.ndarray, rate: int) -> tuple[np.ndarray, np.ndarray]:
    if signal.size < 8:
        return np.zeros(1), np.zeros(1)
    window = np.hanning(signal.size)
    power = np.square(np.abs(np.fft.rfft(window * signal)))
    return np.fft.rfftfreq(signal.size, 1 / rate), power


def _describe(signal: np.ndarray, onset: int, rate: int) -> dict[str, object]:
    regions = []
    for limits in REGIONS:
        audio = _region(signal, onset, rate, limits)
        frequencies, power = _spectrum(audio, rate)
        total = max(float(np.sum(power)), np.finfo(float).tiny)
        fractions = []
        for low, high in BANDS:
            selected = (frequencies >= low) & (frequencies < min(high, rate / 2))
            fractions.append(
                10
                * np.log10(
                    max(float(np.sum(power[selected])), np.finfo(float).tiny) / total
                )
            )
        centroid = float(np.sum(frequencies * power) / total)
        regions.append(
            {
                "seconds": limits,
                "rms_dbfs": 20
                * np.log10(
                    max(float(np.sqrt(np.mean(audio**2))), np.finfo(float).tiny)
                ),
                "centroid_hz": centroid,
                "band_fraction_db": [round(value, 2) for value in fractions],
            }
        )
    return {"onset_samples": onset, "regions": regions}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    arguments = parser.parse_args()
    reference = read_wav(arguments.reference, "mean")
    candidate = resample_audio(
        read_wav(arguments.candidate, "mean"), reference.sample_rate
    )
    reference_onset = detect_impact_onset(reference.samples, reference.sample_rate)
    candidate_onset = detect_impact_onset(candidate.samples, candidate.sample_rate)
    ref_tail = reference.samples[reference_onset:]
    candidate_tail = candidate.samples[candidate_onset:]
    count = min(ref_tail.size, candidate_tail.size)
    denominator = float(np.sum(candidate_tail[:count] ** 2))
    gain = np.sqrt(float(np.sum(ref_tail[:count] ** 2)) / denominator)
    payload = {
        "sample_rate": reference.sample_rate,
        "band_edges_hz": BANDS,
        "candidate_energy_match_gain": gain,
        "reference": _describe(
            reference.samples, reference_onset, reference.sample_rate
        ),
        "candidate": _describe(
            gain * candidate.samples, candidate_onset, candidate.sample_rate
        ),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
