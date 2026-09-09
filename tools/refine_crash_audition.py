"""Bounded crash refinement guided by audition, using only visible DSP controls.

Hold geometry, noisiness, contact and gains fixed. Search six shared edits on
two phase seeds; keep independent seeds for the publication audit. The low
prominence shelf is baked into existing bars, never stored as hidden DSP state.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize
from scipy.stats import qmc

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from crash_beating_common import CrashObjective
from fit_stretched_gong import checkpoint

BOUNDS = {
    "field_phase_bandwidth": (0.002, 0.006),
    "field_packet_spread": (2.4, 3.3),
    "bloom_rate": (2.8, 5.0),
    "body_decay_seconds_0": (18, 28),
    "body_decay_seconds_7": (0.4, 0.85),
    "low_prominence_db": (-7, -1),
}


def audition_parameters(base, position):
    """Map normalized search coordinates to six broad, bounded sound edits."""
    position = np.asarray(position, dtype=float)
    if position.shape != (len(BOUNDS),) or not np.isfinite(position).all():
        raise ValueError("Expected six finite search coordinates")
    if np.any((position < 0) | (position > 1)):
        raise ValueError("Search coordinates must lie in [0, 1]")
    edits = {
        key: float(lo + x * (hi - lo))
        for (key, (lo, hi)), x in zip(BOUNDS.items(), position)
    }
    attenuation = edits.pop("low_prominence_db")
    parameters = dict(base, **edits)
    for i in range(32):
        if base[f"resolved_level_{i}"] <= -71.99:
            continue
        # Full attenuation below 250 Hz, smooth log-frequency fade to 450 Hz.
        f = base[f"resolved_frequency_{i}"]
        weight = np.clip(np.log(450 / f) / np.log(450 / 250), 0, 1)
        parameters[f"resolved_level_{i}"] = max(
            -71, base[f"resolved_level_{i}"] + float(weight) * attenuation
        )
    return parameters


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(reference, renderer.sample_rate, 60)
        args.output.mkdir(parents=True, exist_ok=True)
        seeds = [renderer.metadata["event"]["seed"]]
        seeds.append((seeds[0] + 307) & 0xFFFFFFFF)
        trace, best = [], None

        def evaluate(parameters, position, stage):
            nonlocal best
            rows = [
                objective.components(renderer.render(parameters, 6, s)) for s in seeds
            ]
            scores = [objective.score_components(r) for r in rows]
            row = dict(
                stage=stage,
                position=position,
                score=float(np.mean(scores)),
                seeds=seeds,
                components=rows,
            )
            trace.append(row)
            if best is None or row["score"] < best["score"]:
                best = dict(row, parameters=parameters)
                print(
                    json.dumps({k: v for k, v in row.items() if k != "components"}),
                    flush=True,
                )
            return row["score"]

        evaluate(source["parameters"], None, "unchanged baseline")
        # Separate the user's directional hypothesis from a claim of improvement.
        trials = np.vstack([np.full(6, 0.5), qmc.LatinHypercube(6, seed=29).random(47)])
        scored = []
        for x in trials:
            value = evaluate(
                audition_parameters(source["parameters"], x), x.tolist(), "screen"
            )
            scored.append((value, x))
        initial = min(scored, key=lambda pair: pair[0])[1]
        minimize(
            lambda x: evaluate(
                audition_parameters(source["parameters"], x), x.tolist(), "Powell"
            ),
            initial,
            method="Powell",
            bounds=[(0, 1)] * 6,
            options=dict(maxfev=args.budget, xtol=0.02, ftol=1e-4),
        )
        history = source["history"] + [
            dict(
                stage="audition-directed six-coordinate fit",
                bounds=BOUNDS,
                low_shelf_hz=[250, 450],
                training_seeds=seeds,
                trials=trace,
                accepted=best["position"] is not None,
            )
        ]
        result = checkpoint(
            renderer,
            objective,
            args.output / "candidate",
            "Crash — broader shimmer trial",
            best["parameters"],
            reference,
            history,
        )
        result.seeds = seeds
        result.save()
        verify_candidate(renderer, result.output)
        (args.output / "screen.json").write_text(json.dumps(trace, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=100)
    run(parser.parse_args())
