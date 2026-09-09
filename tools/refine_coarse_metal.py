"""Finish a coarse modal fit with contact and a sparse shared damping curve.

Compare zero/one/two added global T60 knots AFTER the two-endpoint fit. All
painted frequencies and local noisiness stay fixed; prominence has six broad
coordinates. Exact-render Mel plus the existing band-envelope/rise objective
guides scalar fitting; neither is a listening-acceptance test.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint

KNOTS = (120, 400, 1500, 4000, 8000, 15000)


class Score:
    units = "Mel + 0.05 spectral bloom norm"

    def __init__(self, mel, shape):
        self.mel, self.shape = mel, shape
        self.specification = dict(
            mel=mel.specification, shape=shape.specification, shape_weight=0.05
        )

    def score(self, samples):
        return self.mel.score(samples) + 0.05 * np.linalg.norm(
            self.shape.residual(samples)
        )

    def residual(self, samples, regions=None):
        return np.array([np.sqrt(self.score(samples))])


def damping_variant(parameters, frequencies):
    """Build a sparse curve from endpoints, deliberately replacing interior knots."""
    frequencies = np.asarray(frequencies, dtype=float)
    if (
        len(frequencies) > 6
        or not np.isfinite(frequencies).all()
        or np.any(frequencies <= 40)
        or np.any(frequencies >= 15000)
        or np.any(np.diff(frequencies) <= 0)
    ):
        raise ValueError(
            "Use up to six ordered interior frequencies between 40 and 15000 Hz"
        )
    result = dict(parameters)
    erb = lambda f: 21.4 * np.log10(1 + 0.00437 * f)
    for i in range(1, 7):
        result[f"body_decay_active_{i}"] = 0
    for i, frequency in enumerate(frequencies, 1):
        weight = (erb(frequency) - erb(40)) / (erb(15000) - erb(40))
        seconds = np.exp(
            (1 - weight) * np.log(parameters["body_decay_seconds_0"])
            + weight * np.log(parameters["body_decay_seconds_7"])
        )
        result.update(
            {
                f"body_decay_active_{i}": 1,
                f"body_decay_frequency_{i}": frequency,
                f"body_decay_seconds_{i}": float(seconds),
            }
        )
    return result


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
        objective = Score(mel, shape)

        def save(name, parameters, history):
            return checkpoint(
                renderer,
                shape,
                args.output / name,
                f"{args.target.title()} - stretched series, shared decay",
                parameters,
                reference,
                history,
            )

        common = save(
            "contact",
            source["parameters"],
            source["history"] + [dict(parent=str(args.source))],
        )
        if not args.skip_contact:
            common.loss = objective
            refine_scalar(
                common,
                dict(
                    impact_tone_noise=(0.05, 0.95),
                    impact_width=(0.25, 2),
                    impact_noise_tilt=(-6, 18),
                    direct_gain=(0.01, 0.4),
                ),
                budget=180,
                step=0.003,
            )
            common.loss = shape
            polish_coarse(common, KNOTS)
        candidates = []
        for index, frequencies in enumerate(((), (2500,), (500, 2500))):
            if index not in args.variants:
                continue
            search = save(
                f"decay-{index}",
                damping_variant(common.parameters, frequencies),
                common.history
                + [dict(stage="sparse shared curve", frequencies=list(frequencies))],
            )
            bounds = dict(
                bloom_rate=(0.02, 16),
                bloom_energy_acceleration=(0, 1),
                body_brightness=(-36, 6),
                body_excitation_centre=(150, 4500),
                field_turbulence=(0.2, 2.5),
                field_phase_bandwidth=(0.001, 2),
                body_decay_seconds_0=(0.3, 30),
                body_decay_seconds_7=(0.1, 8),
            )
            bounds.update(
                {
                    f"body_decay_seconds_{i}": (0.15, 20)
                    for i in range(1, len(frequencies) + 1)
                }
            )
            scales = {
                key: "log"
                for key in bounds
                if "seconds" in key
                or key
                in ("bloom_rate", "body_excitation_centre", "field_phase_bandwidth")
            }
            search.stage(
                "sparse global damping and bloom",
                bounds,
                28,
                parameter_scales=scales,
                difference_step=0.005,
            )
            polish_coarse(search, KNOTS)
            search.loss = objective
            refine_scalar(
                search, bounds, budget=args.budget, step=0.003, method=args.method
            )
            search.loss = shape
            seed = renderer.metadata["event"]["seed"]
            search.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
            polish_coarse_mel(search, mel, KNOTS)
            value = np.mean(
                [
                    objective.score(search.audio(search.parameters, s))
                    for s in search.seeds
                ]
            )
            print(json.dumps(dict(curve=list(frequencies), score=value)), flush=True)
            candidates.append((value, search))
        _, best = min(candidates, key=lambda row: row[0])
        final = save(
            "candidate", best.parameters, best.history + [dict(parent=str(best.output))]
        )
        verify_candidate(renderer, final.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=260)
    parser.add_argument("--method", choices=["L-BFGS-B", "Powell"], default="Powell")
    parser.add_argument(
        "--skip-contact",
        action="store_true",
        help="Resume from an already contact-refined checkpoint",
    )
    parser.add_argument(
        "--variants", nargs="+", type=int, choices=[0, 1, 2], default=[0, 1, 2]
    )
    run(parser.parse_args())
