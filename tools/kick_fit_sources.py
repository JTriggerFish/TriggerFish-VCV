"""Measure the three unified kick observations without changing their levels."""

import json

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.transforms import StftConfig, stft

ROUTES = ("contact_level", "thump_level", "resonance_level")


def source_audit(renderer, parameters, directory):
    parts, rows = [], []
    directory.mkdir(parents=True, exist_ok=True)
    for route in ROUTES:
        isolated = dict(parameters, **{key: 0 for key in ROUTES if key != route})
        audio = renderer.render(isolated, 1.2)
        parts.append(audio)
        write_wav(directory / f"{route}.wav", AudioBuffer(audio, renderer.sample_rate))
        value = stft(audio, renderer.sample_rate, StftConfig(2048, 256))
        power = value.power
        measurements = []
        for start, end in ((0, 0.03), (0.03, 0.1), (0.1, 0.25)):
            time = (value.times_seconds >= start) & (value.times_seconds < end)
            bands = []
            for low, high in ((0, 120), (120, 1000), (1000, 16000)):
                frequency = (value.frequencies_hz >= low) & (
                    value.frequencies_hz < high
                )
                bands.append(
                    float(
                        10
                        * np.log10(
                            max(power[frequency][:, time].sum(axis=0).mean(), 1e-20)
                        )
                    )
                )
            measurements.append(dict(seconds=[start, end], band_power_db=bands))
        rows.append(dict(route=route, measurements=measurements))
    actual = renderer.render(parameters, 1.2)
    error = float(np.max(np.abs(actual - np.sum(parts, axis=0))))
    report = dict(
        routes=rows,
        bands_hz=[[0, 120], [120, 1000], [1000, 16000]],
        superposition_max_error=error,
    )
    (directory / "sources.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(dict(source_audit=report)), flush=True)
