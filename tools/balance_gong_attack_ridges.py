"""Measured low-ridge level correction; fixed centres, dynamics and master gain.

Use first-400-ms Welch power in +/-6-Hz windows around the measured ridges.
Small bounded observation-bar updates are evaluated in the actual renderer.
This is target-specific voicing, not a general-purpose perceptual objective.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import welch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def ridge_power(samples, rate, centres):
    frequency, power = welch(
        samples[: round(0.4 * rate)], rate, nperseg=16384, nfft=65536
    )
    return np.array(
        [
            10 * np.log10(max(1e-20, power[abs(frequency - f) < 6].sum()))
            for f in centres
        ]
    )


def run(args):
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            SpectralBloomLoss(reference, renderer.sample_rate),
            args.output,
            6,
            "Gong - stretched harmonic starting fit",
        )
        search.parameters = dict(saved["parameters"])
        centres = np.array([121.8, 343.9, 374.1, 553.8])
        indices = [
            min(
                range(32),
                key=lambda i: abs(search.parameters[f"resolved_frequency_{i}"] - f),
            )
            for f in centres
        ]
        if any(
            abs(search.parameters[f"resolved_frequency_{i}"] - f) > 1
            for i, f in zip(indices, centres)
        ):
            raise ValueError("Source must already contain the measured low centres")
        target = ridge_power(reference, renderer.sample_rate, centres)
        search.history.append(
            dict(
                parent=str(args.source.resolve()),
                stage="explicit attack ridge windows",
                centres_hz=centres.tolist(),
                half_width_hz=6,
                attack_seconds=0.4,
                welch_window=16384,
                welch_fft=65536,
            )
        )
        for iteration in range(3):
            before = (
                ridge_power(
                    search.audio(search.parameters), renderer.sample_rate, centres
                )
                - target
            )
            candidate = dict(search.parameters)
            for i, error in zip(indices, before):
                key = f"resolved_level_{i}"
                candidate[key] = float(
                    np.clip(candidate[key] - np.clip(error, -4, 4), -45, 6)
                )
            after = (
                ridge_power(search.audio(candidate), renderer.sample_rate, centres)
                - target
            )
            selected = np.linalg.norm(after) < np.linalg.norm(before)
            search.history.append(
                dict(
                    stage="bounded measured attack-bar correction",
                    iteration=iteration,
                    before_db=before.tolist(),
                    after_db=after.tolist(),
                    selected=bool(selected),
                    trial_parameters=candidate,
                )
            )
            if selected:
                search.parameters = candidate
            print(json.dumps(search.history[-1]), flush=True)
        search.save()
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
