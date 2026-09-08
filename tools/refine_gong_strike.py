"""Refine the gong's initial tonal body without changing its bloom or damping.

Actual Wasm renders; early and full auraloss Mel plus an absolute bloom guard.
Only existing low anchor frequencies, local noisiness and observation levels
move. A harmonic starting hypothesis is compared, never imposed on upper modes.
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


class StrikeLoss(ScalarAudioLoss):
    def __init__(self, reference, rate, baseline):
        self.frames = round(0.4 * rate)
        self.early = AuralossMel(reference[: self.frames], rate)
        self.full = AuralossMel(reference, rate)
        self.bloom = SpectralBloomLoss(reference, rate)
        self.limit = np.linalg.norm(self.bloom.residual(baseline)) + 0.3
        self.specification = dict(
            early=self.early.specification,
            early_seconds=0.4,
            full=self.full.specification,
            weights=[1, 1],
            bloom=self.bloom.specification,
            bloom_limit=float(self.limit),
            penalty="squared excess over limit",
        )

    def score(self, samples):
        excess = max(0.0, np.linalg.norm(self.bloom.residual(samples)) - self.limit)
        return (
            self.early.score(samples[: self.frames])
            + self.full.score(samples)
            + excess**2
        )


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        baseline = renderer.render(saved["parameters"], seconds)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            StrikeLoss(reference, renderer.sample_rate, baseline),
            args.output,
            seconds,
            "Gong - defined low strike",
        )
        search.parameters = saved["parameters"]
        search.history.append(dict(parent=str(args.source.resolve())))
        # The source has ~122,247,375-Hz ridges, but also a strong nonharmonic
        # ~344-Hz ridge. Do not quantize that ridge or the noisy upper spectrum.
        harmonic = dict(search.parameters)
        for index, frequency in [(0, 123.47), (1, 246.94), (4, 370.41)]:
            harmonic[f"resolved_frequency_{index}"] = frequency
        clean = dict(harmonic)
        for index in range(5):
            clean[f"resolved_turbulence_{index}"] = 0.15
        search.screen_candidates(
            "low harmonic hypotheses",
            [("harmonic seed", harmonic), ("clean low harmonics", clean)],
        )
        scalar_stage(
            search,
            "low ridge tuning",
            {
                "resolved_frequency_0": (117, 128),
                "resolved_frequency_1": (239, 255),
                "resolved_frequency_3": (338, 352),
                "resolved_frequency_4": (365, 382),
                "resolved_frequency_6": (535, 560),
            },
            100,
        )
        bounds = {f"resolved_turbulence_{i}": (0, 1) for i in (0, 3, 4)}
        bounds.update(
            {
                f"resolved_level_{i}": (
                    max(-60, search.parameters[f"resolved_level_{i}"] - 6),
                    min(6, search.parameters[f"resolved_level_{i}"] + 10),
                )
                for i in (0, 1, 2, 3, 4, 6)
            }
        )
        scalar_stage(search, "initial ridge balance and clarity", bounds, 180)
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
