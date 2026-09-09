"""Fit the diffusion-only recipe; exact Wasm, fixed gesture and observation level.

Order: shared damping + excitation/transport, observation, packet texture,
joint dynamics, then guarded perceptual observation refinement. Two endpoint
T60s only; fixed frequencies/count, no per-mode decays, no topology switches.
Each solver stage logs bounds and finite-difference influence measurements.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.stats import qmc
import torch

from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_envelope_guard import PerceptualEnvelopeGuard
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from refit_relaxed_metals import checkpoint


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        seconds = 10 if args.target == "crash" else 8
        reference = aligned_reference(renderer, seconds)
        balance = MetallicBalanceLoss(reference, renderer.sample_rate, "erb", True)
        mel = AuralossMel(reference, renderer.sample_rate)
        search = Search(
            renderer,
            balance,
            args.output,
            seconds,
            f"{args.target.title()} — spectral diffusion fit",
            (None,),
        )
        baseline = dict(renderer.initial)
        if args.resume:
            search.parameters = dict(
                verify_candidate(renderer, args.resume)["parameters"]
            )
        else:
            # Visible starting controls, stored in every checkpoint. No hidden
            # recipe coefficients or restored random neighbour exchange.
            search.parameters.update(bloom_energy_acceleration=1, field_wander_hz=0)
            for knot in range(1, 7):
                search.parameters[f"body_decay_active_{knot}"] = 0
        checkpoint(search, "start", baseline, mel)
        damping = dict(body_decay_seconds_0=(0.3, 15), body_decay_seconds_7=(0.15, 10))
        motion = dict(
            damping,
            bloom_rate=(0.01, 16),
            body_brightness=(-60, 12),
            body_excitation_centre=(100, 5000),
        )
        if not args.resume:
            # Screen log-spaced damping/transfer starts before a local Jacobian.
            start = dict(search.parameters)
            candidates = []
            for index, u in enumerate(qmc.LatinHypercube(4, seed=283).random(24)):
                candidates.append(
                    (
                        f"dynamics-{index}",
                        dict(
                            start,
                            body_decay_seconds_0=float(
                                np.exp(np.log(1) + u[0] * np.log(10))
                            ),
                            body_decay_seconds_7=float(
                                np.exp(np.log(0.4) + u[1] * np.log(15))
                            ),
                            bloom_rate=float(np.exp(np.log(0.05) + u[2] * np.log(320))),
                            body_brightness=float(-48 + u[3] * 48),
                        ),
                    )
                )
            search.screen_candidates("diffusion/shared damping starts", candidates)
            search.loss = BandDecayShapeLoss(reference, renderer.sample_rate)
            search.stage("shared damping shape", damping, 16, influence_threshold=0.01)
            search.loss = balance
        search.stage(
            "excitation, diffusion and damping",
            motion,
            20,
            difference_step=0.008,
            influence_threshold=0.01,
        )
        polish_observation_autograd(search, iterations=70)
        checkpoint(search, "dynamics", baseline, mel)
        # Noisiness centre/amount have a multiplicative nullspace. Hold the
        # visible centre fixed; vary amount and slope, never both pivots.
        texture = dict(
            field_turbulence=(0.05, 2),
            field_turbulence_slope=(0, 1),
            field_packet_spread=(0.1, 8),
            field_phase_bandwidth=(0, 1),
        )
        start = dict(search.parameters)
        search.screen_candidates(
            "packet contrast starts",
            [
                (
                    f"noise-{level}-{slope}",
                    dict(
                        start,
                        field_turbulence=level,
                        field_turbulence_slope=slope,
                        field_packet_spread=3,
                        field_phase_bandwidth=0.15,
                    ),
                )
                for level in (0.2, 0.6, 1.2)
                for slope in (0.4, 0.8, 1)
            ],
        )
        search.stage(
            "packet noisiness",
            texture,
            16,
            difference_step=0.015,
            influence_threshold=0.01,
        )
        search.seeds = (None, (renderer.metadata["event"]["seed"] + 21001) & 0xFFFFFFFF)
        search.stage(
            "two-seed dynamics refinement",
            motion,
            16,
            difference_step=0.008,
            influence_threshold=0.01,
        )
        polish_observation_autograd(search, iterations=80)
        checkpoint(search, "balanced", baseline, mel)
        comparator = {
            seed: search.audio(search.parameters, seed).copy() for seed in search.seeds
        }
        search.loss = mel
        guard = lambda basis: PerceptualEnvelopeGuard(
            basis, reference, renderer.sample_rate, comparator, tolerance_db=0.5
        )
        polish_observation_autograd(search, iterations=60, constraint_factory=guard)
        search.loss = balance
        checkpoint(search, "perceptual", baseline, mel)
        print(
            json.dumps(dict(done=args.target, evaluations=search.evaluations)),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", type=Path)
    run(parser.parse_args())
