"""Library Mel refinement of four broad observation coordinates, with bloom guard."""

import argparse
import os
from pathlib import Path

import torch

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel as polish
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        search = checkpoint(
            renderer,
            shape,
            args.output,
            "Gong protected harmonic core",
            saved["parameters"],
            reference,
            saved["history"] + [dict(parent=str(args.source))],
        )
        search.seeds = (1675, 1776)
        polish(search, AuralossMel(reference, renderer.sample_rate))
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
