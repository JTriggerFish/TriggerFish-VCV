"""Use auraloss mel while guarding the already fitted onset/rise/decay cells."""

import argparse
import os
from pathlib import Path

import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomGuard
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.source)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            AuralossMel(reference, renderer.sample_rate),
            args.output,
            seconds,
            args.target + " - bloom and spectral refinement",
            seeds=tuple(args.seeds) if args.seeds else (None,),
        )
        search.parameters = saved["parameters"]
        search.history.append(
            dict(parent=str(args.source.resolve()), guard=shape.specification)
        )
        fixed = [
            f"resolved_level_{i}"
            for i in range(32)
            if args.fixed_below
            and search.parameters[f"resolved_frequency_{i}"] < args.fixed_below
        ]
        polish_observation_autograd(
            search,
            iterations=35,
            constraint_factory=lambda basis: SpectralBloomGuard(basis, shape, 1.0),
            fixed_keys=fixed,
        )
        search.save()
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+")
    parser.add_argument(
        "--fixed-below",
        type=float,
        default=0,
        help="Hold measured observation bars below this frequency fixed",
    )
    run(parser.parse_args())
