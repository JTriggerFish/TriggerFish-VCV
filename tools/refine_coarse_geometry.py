"""Test measured low/mid anchor corrections around a fixed stretched series.

No mode-frequency optimizer: root and optional index=Hz corrections are
explicit inputs. Corrected centres may deviate at most 10% from the scaled
series. All prominence fitting retains the same six broad coordinates.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.coarse_observation_fit import (
    interpolation_weights,
    polish_coarse,
)
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_coarse_metal import KNOTS


def geometry(parameters, root, corrections):
    result = dict(parameters)
    indices = [i for i in range(32) if parameters[f"resolved_level_{i}"] > -71.99]
    old = np.array([parameters[f"resolved_frequency_{i}"] for i in indices])
    levels = np.array([parameters[f"resolved_level_{i}"] for i in indices])
    coefficients = np.linalg.lstsq(
        interpolation_weights(old, KNOTS), 10 ** (levels / 20), rcond=None
    )[0]
    frequencies = old * root / old[0]
    changes = []
    for index, frequency in corrections.items():
        slot = indices.index(index)
        if not np.isfinite(frequency) or abs(frequency / frequencies[slot] - 1) > 0.1:
            raise ValueError("Measured anchor correction exceeds 10% of scaled series")
        changes.append(
            dict(index=index, series_hz=float(frequencies[slot]), measured_hz=frequency)
        )
        frequencies[slot] = frequency
    if (
        not np.isfinite(frequencies).all()
        or np.any(np.diff(frequencies) <= 0)
        or not 1 <= frequencies[0] < frequencies[-1] <= 15000
    ):
        raise ValueError("Corrected series must remain ordered inside the UI range")
    amplitudes = interpolation_weights(frequencies, KNOTS) @ coefficients
    if np.any(amplitudes <= 0):
        raise ValueError("Source does not support a positive broad observation curve")
    for index, frequency, amplitude in zip(indices, frequencies, amplitudes):
        result[f"resolved_frequency_{index}"] = float(frequency)
        result[f"resolved_level_{index}"] = float(20 * np.log10(amplitude))
    return result, changes


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.source)
        corrections = {
            int(pair.split("=")[0]): float(pair.split("=")[1]) for pair in args.anchor
        }
        parameters, changes = geometry(saved["parameters"], args.root, corrections)
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        mel = AuralossMel(reference, renderer.sample_rate)
        search = checkpoint(
            renderer,
            shape,
            args.output,
            f"{args.target.title()} - stretched series with measured low-mid anchors",
            parameters,
            reference,
            saved["history"]
            + [
                dict(
                    stage="explicit low-mid geometry correction",
                    parent=str(args.source),
                    root_hz=args.root,
                    anchors=changes,
                )
            ],
        )
        polish_coarse(search, KNOTS)
        seed = renderer.metadata["event"]["seed"]
        search.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
        polish_coarse_mel(search, mel, KNOTS)
        verify_candidate(renderer, args.output)
        audio = search.audio(search.parameters)
        print(
            json.dumps(
                dict(
                    mel=mel.score(audio),
                    shape=float(np.linalg.norm(shape.residual(audio))),
                )
            ),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--root", type=float, required=True)
    parser.add_argument(
        "--anchor", action="append", default=[], help="Zero-based index=measured Hz"
    )
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
