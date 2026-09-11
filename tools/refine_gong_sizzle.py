"""Screen packet texture against regional brightness/modulation on exact WASM.

Keep every centre, bar, damping point, excitation and gain fixed in this pass.
Report per-seed diagnostics; do not accept a mean score as listening approval.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.upper_sizzle_loss import UpperSizzleLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def candidates(base):
    yield "baseline", dict(base)
    for key, values in {
        "field_packet_spread": [1, 1.8, 4, 6],
        "field_satellite_density": [0.15, 0.3, 0.7, 1],
        "field_phase_bandwidth": [0, 0.01, 0.07, 0.15],
        "field_phase_tilt": [-1, 0, 0.5],
        "field_turbulence_slope": [0.2, 0.6, 0.85],
    }.items():
        for value in values:
            yield f"{key}={value}", dict(base, **{key: value})
    # Smooth upper-only allocation/profile edits leave the low handles alone.
    for allocation in (0.25, 2, 4):
        p = dict(base)
        for i in range(32):
            if p[f"resolved_frequency_{i}"] >= 3000:
                p[f"resolved_allocation_{i}"] = allocation
        yield f"upper allocation {allocation}", p


def joint_candidates(base, stage):
    """Coordinate grids around the incumbent, not unconstrained per-mode fitting."""
    grids = [
        {
            "field_phase_bandwidth": [0.035, 0.1, 0.25, 0.5],
            "field_phase_tilt": [0.5, 1, 1.5],
        },
        {
            "field_packet_spread": [1.8, 2.7, 4, 6],
            "field_satellite_density": [0.45, 0.7, 1],
        },
        {"field_turbulence_slope": [0.2, 0.4, 0.7, 1], "bloom_rate": [4, 6, 8, 12]},
    ]
    from itertools import product

    grid = grids[stage]
    for values in product(*grid.values()):
        change = dict(zip(grid, values))
        yield str(change), dict(base, **change)


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        baseline_path = args.output / "baseline.json"
        if baseline_path.exists():
            base = json.loads(baseline_path.read_text())
        else:
            base = dict(r.initial)
            baseline_path.write_text(json.dumps(base, indent=2))
        ref = aligned_reference(r, 6)
        loss = UpperSizzleLoss(ref, r.sample_rate)
        seeds = (1675, 1982, 2586)
        originals = {s: loss.measure(r.render(base, 6, s)) for s in seeds}
        rows = []

        def evaluate(name, p):
            metrics = [loss.diagnostics(r.render(p, 6, s), originals[s]) for s in seeds]
            scores = [m["score"] for m in metrics]
            row = dict(
                name=name,
                parameters=p,
                metrics=metrics,
                score=float(np.mean(scores) + 0.5 * max(scores)),
            )
            rows.append(row)
            # Hard guard protects every measured low-body region and seed.
            row["eligible"] = max(m["low_max_change_db"] for m in metrics) <= 2
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            print(
                json.dumps(
                    {k: v for k, v in row.items() if k not in ("parameters", "metrics")}
                ),
                flush=True,
            )

        if args.joint:
            evaluate("baseline", base)
            for stage in range(3):
                incumbent = min(
                    (x for x in rows if x["eligible"]), key=lambda x: x["score"]
                )
                for name, p in joint_candidates(incumbent["parameters"], stage):
                    evaluate(name, p)
        else:
            for name, p in candidates(base):
                evaluate(name, p)
        selected = min((x for x in rows if x["eligible"]), key=lambda x: x["score"])
        search = checkpoint(
            r,
            loss,
            args.output / "candidate",
            "Gong — brighter sizzle",
            selected["parameters"],
            ref,
            [
                dict(
                    stage="texture screen",
                    selected=selected["name"],
                    objective=loss.specification,
                    seeds=seeds,
                    baseline=rows[0]["score"],
                    score=selected["score"],
                )
            ],
        )
        verify_candidate(r, search.output)
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--joint", action="store_true")
    run(parser.parse_args())
