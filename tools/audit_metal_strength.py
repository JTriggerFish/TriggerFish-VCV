"""Exact-workbench, fixed-gain strength and early-bloom ablations.

STFT power is integrated in physical frequency bands; no per-render gain or
colour normalization. This is diagnostic evidence, not an acceptance loss.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import stft

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

BANDS = ((40, 300), (300, 1000), (1000, 3000), (3000, 6000), (6000, 16000))
REGIONS = ((0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1), (1, 2), (2, 4))


def measure(samples, rate):
    frequencies, times, spectrum = stft(
        samples,
        rate,
        nperseg=4096,
        noverlap=4096 - round(0.01 * rate),
        boundary="zeros",
    )
    power = np.array(
        [
            np.sum(abs(spectrum[(frequencies >= lo) & (frequencies < hi)]) ** 2, axis=0)
            for lo, hi in BANDS
        ]
    )
    windows = np.array(
        [np.mean(power[:, (times >= lo) & (times < hi)], axis=1) for lo, hi in REGIONS]
    ).T
    return dict(
        peak_db=float(20 * np.log10(max(1e-20, abs(samples).max()))),
        band_region_db=(10 * np.log10(np.maximum(windows, 1e-20))).tolist(),
        band_peak_seconds=times[np.argmax(power, axis=1)].tolist(),
    )


def run(args):
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        parameters = renderer.initial
        if args.parameters:
            saved = json.loads(args.parameters.read_text())
            parameters = saved.get("parameters", saved)
        output = args.output
        output.mkdir(parents=True, exist_ok=True)
        reference = aligned_reference(renderer, 6)
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        standard = renderer.render(parameters, 6)
        write_wav(output / "candidate.wav", AudioBuffer(standard, renderer.sample_rate))
        rows = {
            "reference": measure(reference, renderer.sample_rate),
            "standard": measure(standard, renderer.sample_rate),
        }
        ablations = {
            "current": {},
            "no-transfer": {"bloom_rate": 0},
            "no-phase-noise": {"field_phase_bandwidth": 0},
            "no-transfer-or-phase": {"bloom_rate": 0, "field_phase_bandwidth": 0},
        }
        for name, override in ablations.items():
            for strength in sorted(
                set((0.25, 0.5, 0.75, 1, renderer.metadata["event"]["strength"]))
            ):
                key = f"{name}-{strength:.3f}"
                samples = renderer.decode(
                    renderer.request(
                        command="renderSequence",
                        parameters=dict(parameters, **override),
                        seconds=6,
                        hits=[dict(time=0, strength=strength)],
                    )["pcm"]
                )
                rows[key] = measure(samples, renderer.sample_rate)
                write_wav(
                    output / (key + ".wav"), AudioBuffer(samples, renderer.sample_rate)
                )
                print(json.dumps({key: rows[key]}), flush=True)
        report = dict(
            renderer=renderer.metadata,
            parameters=parameters,
            bands=BANDS,
            regions=REGIONS,
            rows=rows,
        )
        (output / "audit.json").write_text(json.dumps(report, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parameters", type=Path)
    run(parser.parse_args())
