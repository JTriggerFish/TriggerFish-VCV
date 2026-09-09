"""Compare global packet textures with fixed geometry and broad prominence only."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_coarse_metal import KNOTS


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        mel = AuralossMel(reference, renderer.sample_rate)

        def score(search):
            audio = search.audio(search.parameters)
            return mel.score(audio) + 0.05 * np.linalg.norm(shape.residual(audio))

        candidates = []
        for level in args.levels:
            for spread in args.spreads:
                parameters = dict(
                    source["parameters"],
                    field_turbulence=level,
                    field_turbulence_slope=1,
                    field_packet_spread=spread,
                    field_phase_bandwidth=args.bandwidth,
                )
                settings = dict(
                    level=level, spread=spread, bandwidth=args.bandwidth, slope=1
                )
                search = checkpoint(
                    renderer,
                    shape,
                    args.output / f"texture-{level:g}-{spread:g}",
                    f"{args.target.title()} - stretched series, global texture fit",
                    parameters,
                    reference,
                    source["history"]
                    + [
                        dict(
                            stage="global texture screen",
                            settings=settings,
                            parent=str(args.source),
                        )
                    ],
                )
                polish_coarse(search, KNOTS)
                value = score(search)
                candidates.append((value, search))
                print(json.dumps(dict(**settings, score=float(value))), flush=True)
        finished = []
        for _, search in sorted(candidates, key=lambda item: item[0])[:2]:
            seed = renderer.metadata["event"]["seed"]
            search.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
            polish_coarse_mel(search, mel, KNOTS)
            finished.append((score(search), search))
        _, best = min(finished, key=lambda item: item[0])
        final = checkpoint(
            renderer,
            shape,
            args.output / "candidate",
            best.name,
            best.parameters,
            reference,
            best.history + [dict(parent=str(best.output))],
        )
        verify_candidate(renderer, final.output)
        print(json.dumps(dict(score=float(score(final)))), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=float, nargs="+", default=[0.6, 0.9, 1.2])
    parser.add_argument("--spreads", type=float, nargs="+", default=[0.35, 0.7, 1.4])
    parser.add_argument("--bandwidth", type=float, default=0.3)
    run(parser.parse_args())
