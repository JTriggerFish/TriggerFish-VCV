"""Compare broad prominence-curve centres without changing mode frequencies."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear
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


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.source)
        original = saved["parameters"]
        indices = [i for i in range(32) if original[f"resolved_level_{i}"] > -71.99]
        frequencies = [original[f"resolved_frequency_{i}"] for i in indices]
        amplitudes = 10 ** (
            np.array([original[f"resolved_level_{i}"] for i in indices]) / 20
        )
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        mel = AuralossMel(reference, renderer.sample_rate)
        best = None
        for centre in args.centres:
            knots = (120, 400, centre, 4000, 8000, 15000)
            weights = interpolation_weights(frequencies, knots)
            projection = lsq_linear(
                weights, amplitudes, bounds=(10 ** (-45 / 20), 10 ** (6 / 20))
            ).x
            parameters = dict(original)
            for index, amplitude in zip(indices, weights @ projection):
                parameters[f"resolved_level_{index}"] = float(20 * np.log10(amplitude))
            search = checkpoint(
                renderer,
                shape,
                args.output / f"centre-{centre:g}",
                f"{args.target.title()} - stretched series, broad prominence fit",
                parameters,
                reference,
                saved["history"]
                + [
                    dict(
                        stage="explicit broad curve projection",
                        parent=str(args.source),
                        knots_hz=knots,
                    )
                ],
            )
            polish_coarse(search, knots)
            seed = renderer.metadata["event"]["seed"]
            search.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
            polish_coarse_mel(search, mel, knots)
            scores = [
                (
                    mel.score(search.audio(search.parameters, s)),
                    float(
                        np.linalg.norm(
                            shape.residual(search.audio(search.parameters, s))
                        )
                    ),
                )
                for s in search.seeds
            ]
            score = float(np.mean([a + 0.05 * b for a, b in scores]))
            print(
                json.dumps(dict(centre=centre, scores=scores, score=score)), flush=True
            )
            if best is None or score < best[0]:
                best = (score, search)
        search = best[1]
        final = checkpoint(
            renderer,
            shape,
            args.output / "candidate",
            search.name,
            search.parameters,
            reference,
            search.history + [dict(parent=str(search.output))],
        )
        verify_candidate(renderer, final.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--centres", nargs="+", type=float, default=[1800, 2100, 2400])
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
