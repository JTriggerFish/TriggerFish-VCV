"""Small visible-frequency/texture refinement, constrained by bloom accuracy.

Use exact Wasm with bounded Powell steps. No extra modes, damping multipliers,
output gain, velocity curve or hidden corrective filter is fitted.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel, ScalarAudioLoss
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from refine_metal_perceptual import scalar_stage


class GuardedMel(ScalarAudioLoss):
    def __init__(self, reference, rate, baseline):
        self.mel = AuralossMel(reference, rate)
        self.bloom = SpectralBloomLoss(reference, rate)
        self.limit = np.linalg.norm(self.bloom.residual(baseline)) + 0.25
        self.specification = dict(
            mel=self.mel.specification,
            bloom=self.bloom.specification,
            bloom_rms_limit=float(self.limit),
            penalty="squared excess over limit",
        )

    def score(self, samples):
        excess = max(0.0, np.linalg.norm(self.bloom.residual(samples)) - self.limit)
        return self.mel.score(samples) + excess**2


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.source)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        baseline = renderer.render(saved["parameters"], seconds)
        loss = GuardedMel(reference, renderer.sample_rate, baseline)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer, loss, args.output, seconds, args.target + " - ridge refinement"
        )
        search.parameters = saved["parameters"]
        search.history.append(dict(parent=str(args.source.resolve())))
        # Only the six most prominent low/mid receiving handles. High-band
        # noise is not a set of tonal ridge frequencies to align independently.
        active = [
            i
            for i in range(32)
            if search.parameters[f"resolved_level_{i}"] > -20
            and search.parameters[f"resolved_frequency_{i}"] < 5000
        ]
        active.sort(
            key=lambda i: search.parameters[f"resolved_level_{i}"], reverse=True
        )
        bounds = {
            f"resolved_frequency_{i}": (
                max(40, search.parameters[f"resolved_frequency_{i}"] * 0.92),
                min(15000, search.parameters[f"resolved_frequency_{i}"] * 1.08),
            )
            for i in active[:6]
        }
        # The exposed 450-Hz handle is the conspicuous low ridge in this crash.
        if args.target == "crash":
            bounds["resolved_frequency_24"] = (390, 470)
        scalar_stage(search, "resolved ridge placement", bounds, 160)
        scalar_stage(
            search,
            "packet spectral width",
            dict(
                field_turbulence=(0.1, 2),
                field_packet_spread=(0.2, 5),
                field_phase_bandwidth=(0.001, 1.2),
            ),
            90,
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
    run(parser.parse_args())
