"""Refine shared gong dynamics after an explicitly reviewed modal start.

No mode placement/level, output gain, EQ, gesture, allocation or extra decay
knots are fitted. Two fixed seeds train eight shared coordinates with bounded
Powell searches; local finite differences and every accepted stage are logged.
The scalar ranking combines fixed-reference Mel and band envelope errors.
It is a search criterion, not perceptual approval: inspect STFTs and audit
additional seeds before publishing to the main workbench.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


class GongEnvelopeLoss:
    """Keep absolute spectral level and relative decay as separate terms."""

    units = "Mel + 0.04 * bloom/rise dB norm + 0.025 * decay dB error"

    def __init__(self, reference, rate):
        self.mel = ReferenceFloorMel(reference, rate, 60)
        self.bloom = SpectralBloomLoss(reference, rate)
        self.decay = BandDecayShapeLoss(reference, rate)
        self.specification = dict(
            mel=self.mel.specification,
            bloom=self.bloom.specification,
            decay=self.decay.specification,
            weights=dict(mel=1, bloom=0.04, decay=0.025),
            normalization=False,
        )

    def score(self, samples):
        return (
            self.mel.score(samples)
            + 0.04 * np.linalg.norm(self.bloom.residual(samples))
            + 0.025 * self.decay.diagnostics(samples)["shape_error_db"]
        )


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        loss = GongEnvelopeLoss(reference, renderer.sample_rate)
        search = checkpoint(
            renderer,
            loss,
            args.output,
            "Gong — low ringing and shared decay",
            saved["parameters"],
            reference,
            saved["history"] + [dict(parent=str(args.source))],
        )
        seed = renderer.metadata["event"]["seed"]
        search.seeds = (seed, seed + 307)
        bounds = dict(
            bloom_rate=(1, 8),
            bloom_energy_acceleration=(0, 0.3),
            bloom_energy_sensitivity=(0, 2),
            body_brightness=(-48, -12),
            field_packet_spread=(2, 5),
            body_decay_seconds_0=(5, 14),
            body_decay_seconds_7=(0.8, 4),
            field_wander_hz=(0.1, 8),
        )
        refine_scalar(search, bounds, budget=args.budget, step=0.01, method="Powell")
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=200)
    run(parser.parse_args())
