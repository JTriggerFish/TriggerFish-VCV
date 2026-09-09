"""Broad observation and optional single-knot damping refinement with texture checks."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import lsq_linear
from triggerfish_percussion.coarse_observation_fit import (
    interpolation_weights,
    polish_coarse,
)
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from fit_structured_metal_texture import Objective
from refine_coarse_metal import damping_variant


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        objective = Objective(reference, renderer.sample_rate)
        seed = renderer.metadata["event"]["seed"]
        seeds = (seed, (seed + 101) & 0xFFFFFFFF)

        def score(p):
            return float(
                np.mean([objective.score(renderer.render(p, 6, s)) for s in seeds])
            )

        original = source["parameters"]
        best, best_score = original, score(original)
        best_history = source["history"]
        active = [i for i in range(32) if original[f"resolved_level_{i}"] > -71.99]
        f = [original[f"resolved_frequency_{i}"] for i in active]
        a = 10 ** (np.array([original[f"resolved_level_{i}"] for i in active]) / 20)
        knots = (120, 400, 1000, 2500, 6500, 15000)
        weights = interpolation_weights(f, knots)
        projection = lsq_linear(weights, a, bounds=(10 ** (-45 / 20), 10 ** (6 / 20))).x
        projected = dict(original)
        for i, level in zip(active, 20 * np.log10(weights @ projection)):
            projected[f"resolved_level_{i}"] = float(level)
        for interior in ([], [3000]):
            p = damping_variant(projected, interior)
            trial = checkpoint(
                renderer,
                objective.shape,
                args.output / f"knots-{len(interior)}",
                args.target.title() + " — structured beating fit",
                p,
                reference,
                source["history"]
                + [
                    dict(
                        stage="six broad prominence coordinates",
                        knots_hz=knots,
                        interior_t60_hz=interior,
                    )
                ],
            )
            trial.seeds = seeds
            polish_coarse(trial, knots)
            trial.loss = objective
            bounds = dict(
                body_decay_seconds_0=(0.3, 30),
                body_decay_seconds_7=(0.1, 15),
                bloom_rate=(0.1, 16),
                bloom_energy_acceleration=(0, 0.25),
            )
            if interior:
                bounds["body_decay_seconds_1"] = (0.1, 20)
            refine_scalar(
                trial, bounds, budget=args.budget, step=0.004, method="Powell"
            )
            before = dict(trial.parameters)
            before_score = score(before)
            trial.loss = objective.shape
            polish_coarse(trial, knots)
            trial.loss = objective
            after_score = score(trial.parameters)
            if after_score > before_score:
                trial.parameters = before
            value = min(after_score, before_score)
            trial.history.append(
                dict(
                    stage="texture guarded prominence",
                    before=before_score,
                    proposal=after_score,
                    accepted=after_score <= before_score,
                )
            )
            trial.save()
            print(json.dumps(dict(interior=interior, score=value)), flush=True)
            # Require a meaningful improvement before adding the extra knot.
            if value < best_score - (0.015 if interior else 0):
                best, best_score, best_history = trial.parameters, value, trial.history
        result = checkpoint(
            renderer,
            objective,
            args.output / "candidate",
            args.target.title() + " — structured beating fit",
            best,
            reference,
            best_history
            + [dict(stage="accepted sparse refinement", score=best_score, seeds=seeds)],
        )
        verify_candidate(renderer, result.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=140)
    run(parser.parse_args())
