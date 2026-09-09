"""Test an explicitly fast, energy-dependent crash front with low fixed blur.

Keeps the modal family and packet texture fixed. This is a bounded hypothesis,
not a new DSP path and not automatic preset publication.
"""

import argparse
import json
import os
from pathlib import Path

import torch

from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_crash_beating import CrashObjective, KNOTS, guarded_polish


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(
            reference, renderer.sample_rate, args.reference_floor_db
        )
        p = dict(
            source["parameters"],
            bloom_rate=8,
            bloom_energy_acceleration=0.2,
            body_brightness=-12,
            body_excitation_centre=800,
            field_phase_bandwidth=0.0005,
        )
        search = checkpoint(
            renderer,
            objective,
            args.output,
            "Crash — fast front trial",
            p,
            reference,
            source["history"]
            + [
                dict(
                    stage="fast front hypothesis",
                    fixed_phase_blur=0.0005,
                    initial_overrides={
                        k: p[k]
                        for k in (
                            "bloom_rate",
                            "bloom_energy_acceleration",
                            "body_brightness",
                            "body_excitation_centre",
                        )
                    },
                )
            ],
        )
        search.loss = objective.shape
        polish_coarse(search, KNOTS)
        search.loss = objective
        refine_scalar(
            search,
            dict(
                bloom_rate=(3, 16),
                bloom_energy_acceleration=(0.1, 0.5),
                body_brightness=(-36, 0),
                body_excitation_centre=(300, 4000),
                body_decay_seconds_0=(0.3, 30),
                body_decay_seconds_7=(0.1, 10),
            ),
            budget=args.budget,
            step=0.004,
            method="Powell",
        )
        guarded_polish(search, objective)
        verify_candidate(renderer, search.output)
        print(
            json.dumps(objective.components(search.audio(search.parameters))),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=220)
    parser.add_argument("--reference-floor-db", type=float)
    run(parser.parse_args())
