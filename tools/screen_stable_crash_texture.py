"""Test smooth drift and broad allocation tilt without stochastic phase blur.

Uses existing controls only. A whole smooth allocation curve replaces manual
per-handle fitting; drift is a constructive alternative, not a physical law.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_crash_beating import CrashObjective, KNOTS, guarded_polish


def allocation_tilt(parameters, exponent):
    if not np.isfinite(exponent) or not 0 <= exponent <= 1:
        raise ValueError("Allocation tilt must be in [0, 1]")
    active = [i for i in range(32) if parameters[f"resolved_level_{i}"] > -71.99]
    result = dict(parameters)
    if not active:
        return result
    maximum = max(parameters[f"resolved_frequency_{i}"] for i in active)
    for i in active:
        result[f"resolved_allocation_{i}"] = (
            1
            if exponent == 0
            else max(
                0.1, 4 * (parameters[f"resolved_frequency_{i}"] / maximum) ** exponent
            )
        )
    return result


def screen(renderer, base, objective):
    rows = []
    for tilt in (0, 0.5, 1):
        for depth, rate in (
            (0, 2),
            (0.03, 2),
            (0.03, 10),
            (0.1, 2),
            (0.1, 10),
            (0.3, 2),
            (0.3, 10),
        ):
            p = allocation_tilt(base, tilt)
            p.update(
                field_phase_bandwidth=0, field_wander_hz=depth, field_wander_rate=rate
            )
            audio = renderer.render(p, 6)
            row = dict(
                tilt=tilt,
                depth=depth,
                rate=rate,
                parameters=p,
                score=objective.score(audio),
                components=objective.components(audio),
            )
            rows.append(row)
            print(
                json.dumps({k: v for k, v in row.items() if k != "parameters"}),
                flush=True,
            )
    return sorted(rows, key=lambda row: row["score"])


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(reference, renderer.sample_rate)
        rows = screen(renderer, source["parameters"], objective)
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
        baseline = checkpoint(
            renderer,
            objective,
            args.output / "baseline",
            "Crash — stable starting point",
            source["parameters"],
            reference,
            source["history"],
        )
        best = objective.score(baseline.audio(baseline.parameters)), baseline
        for index, row in enumerate(rows[:3]):
            trial = checkpoint(
                renderer,
                objective,
                args.output / f"trial-{index}",
                "Crash — stable texture trial",
                row["parameters"],
                reference,
                source["history"] + [dict(stage="stable texture screen", **row)],
            )
            trial.loss = objective.shape
            polish_coarse(trial, KNOTS)
            trial.loss = objective
            guarded_polish(trial, objective)
            value = objective.score(trial.audio(trial.parameters))
            if value < best[0]:
                best = value, trial
        _, chosen = best
        result = checkpoint(
            renderer,
            objective,
            args.output / "candidate",
            "Crash — stable texture trial",
            chosen.parameters,
            reference,
            chosen.history,
        )
        verify_candidate(renderer, result.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
