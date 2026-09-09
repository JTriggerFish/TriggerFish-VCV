"""Rebalance a metallic fit across seeds, then constrain perceptual polishing.

This uses the existing absolute-energy and relative-decay guards rather than
allowing a perceptual objective to trade away body energy unnoticed.
"""

import argparse
import json
import os
from pathlib import Path

import torch

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_envelope_guard import PerceptualEnvelopeGuard
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from refine_metal_perceptual import checkpoint


def fit(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.resume)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        seed = renderer.metadata["event"]["seed"]
        seeds = (None, (seed + args.seed_offset) & 0xFFFFFFFF)
        args.output.mkdir(parents=True, exist_ok=True)
        balance = MetallicBalanceLoss(reference, renderer.sample_rate, "erb", True)
        mel = AuralossMel(reference, renderer.sample_rate)
        search = Search(
            renderer,
            balance,
            args.output,
            seconds,
            f"{args.target.title()} — relaxed reference refinement",
            seeds,
        )
        search.parameters = dict(saved["parameters"])
        if args.cascade_min:
            search.parameters["bloom_rate"] = max(
                args.cascade_min, search.parameters["bloom_rate"]
            )
        baseline = dict(search.parameters)
        search.history.append(
            dict(parent=str(args.resume.resolve()), training_seeds=seeds)
        )
        polish_observation_autograd(search, iterations=100)
        decay = {"body_decay_seconds_0": (0.2, 30), "body_decay_seconds_7": (0.15, 8)}
        for knot in range(1, 7):
            if search.parameters[f"body_decay_active_{knot}"] >= 0.5:
                decay[f"body_decay_seconds_{knot}"] = (0.2, 25)
        transport = dict(decay, bloom_rate=(max(0.01, args.cascade_min), 8))
        if not args.cascade_min:
            transport["bloom_energy_acceleration"] = (0, 1)
        search.stage(
            "multi-seed shared damping and transport",
            transport,
            25,
            difference_step=0.005,
            influence_threshold=0.01,
        )
        polish_observation_autograd(search, iterations=100)
        # Keep a complete balance candidate before changing the objective.
        search.loss = mel
        checkpoint(search, "balanced", baseline, reference)
        comparator = {
            seed: search.audio(search.parameters, seed).copy() for seed in seeds
        }
        factory = lambda basis: PerceptualEnvelopeGuard(
            basis, reference, renderer.sample_rate, comparator, tolerance_db=0.5
        )
        polish_observation_autograd(search, iterations=100, constraint_factory=factory)
        checkpoint(search, "guarded", baseline, reference)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed-offset", type=int, required=True)
    parser.add_argument("--cascade-min", type=float, default=0)
    args = parser.parse_args()
    if not 0 < args.seed_offset < 2**32:
        parser.error("Expected a nonzero 32-bit training seed offset")
    if not 0 <= args.cascade_min < 8:
        parser.error("Cascade minimum must be in [0, 8)")
    fit(args)
