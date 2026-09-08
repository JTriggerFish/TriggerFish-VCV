"""Refine rise/peak/decay using the existing tested regional-energy objective.

No new DSP controls or ad-hoc onset delays. Explicit equal-weight time cells
prevent a missing high-band rise being lost in a long whole-tail score.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.stats import qmc
import torch

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.regional_energy_loss import RegionalEnergyLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from refit_relaxed_metals import checkpoint
from fit_regional_observation import polish


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.source)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        regions = (
            (0, 0.1),
            (0.1, 0.3),
            (0.3, 0.6),
            (0.6, 1),
            (1, 1.5),
            (1.5, 2.5),
            (2.5, 4),
            (4, 6),
        )
        bands = (
            (100, 300),
            (300, 700),
            (700, 1500),
            (1500, 3000),
            (3000, 6000),
            (6000, 16000),
        )
        loss = RegionalEnergyLoss(reference, renderer.sample_rate, bands, regions)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            loss,
            args.output,
            seconds,
            args.target.title() + " — diffusion bloom refinement",
            (None,),
        )
        search.parameters = dict(saved["parameters"])
        baseline = dict(search.parameters)
        search.history.append(
            dict(parent=str(args.source.resolve()), measurement=loss.specification)
        )
        bounds = dict(
            bloom_rate=(0.005, 16),
            body_brightness=(-60, 0),
            body_excitation_centre=(100, 4000),
            field_turbulence=(0.05, 2),
            field_turbulence_slope=(0, 1),
            field_phase_bandwidth=(0, 1),
            body_decay_seconds_0=(0.5, 12),
            body_decay_seconds_7=(0.15, 8),
        )
        if args.clean_source:
            # Diagnostic constrained branch: test whether broadband phase
            # skirts and direct high excitation conceal the transport bloom.
            # These remain ordinary visible values, not DSP overrides.
            search.parameters.update(
                field_phase_bandwidth=0, body_brightness=-36, body_excitation_centre=300
            )
            bounds.update(body_brightness=(-72, -18), body_excitation_centre=(40, 1000))
            del bounds["field_phase_bandwidth"]
        # Deliberately cross the steep/clean excitation alternatives rather
        # than relying solely on a Jacobian around the old immediate wash.
        candidates = []
        for i, u in enumerate(qmc.LatinHypercube(5, seed=1773).random(32)):
            candidates.append(
                (
                    f"rise-{i}",
                    dict(
                        search.parameters,
                        bloom_rate=float(0.005 * 3200 ** u[0]),
                        body_brightness=float(-60 * u[1]),
                        body_excitation_centre=float(100 * 40 ** u[2]),
                        field_turbulence=float(0.1 + 1.2 * u[3]),
                        field_turbulence_slope=float(0.4 + 0.6 * u[4]),
                    ),
                )
            )
            if args.clean_source:
                candidates[-1][1].update(
                    body_brightness=float(-72 + 54 * u[1]),
                    body_excitation_centre=float(40 * 25 ** u[2]),
                )
        search.screen_candidates("rise/decay alternatives", candidates)
        search.stage(
            "regional rise and decay",
            bounds,
            22,
            difference_step=0.01,
            influence_threshold=0.01,
        )
        polish(search, 120)
        search.seeds = (None, (renderer.metadata["event"]["seed"] + 21001) & 0xFFFFFFFF)
        search.stage(
            "regional two-seed refinement",
            bounds,
            14,
            difference_step=0.008,
            influence_threshold=0.01,
        )
        polish(search, 120)
        regional = loss.diagnostics(search.audio(search.parameters))
        search.history.append(dict(stage="final regional review", diagnostics=regional))
        search.loss = MetallicBalanceLoss(reference, renderer.sample_rate, "erb", True)
        checkpoint(
            search, "candidate", baseline, AuralossMel(reference, renderer.sample_rate)
        )
        print(json.dumps(dict(done=args.target, regional=regional)), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--clean-source", action="store_true")
    run(parser.parse_args())
