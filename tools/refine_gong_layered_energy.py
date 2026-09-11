"""Check whether upper observation saturation needs more transported energy.

Screen only existing excitation/transport controls; then refit a broad painted
balance. No output gain adjustment or additional DSP control is introduced.
"""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from refine_gong_layered_timing import polish


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = json.loads(args.source.read_text())["parameters"]
        ref = aligned_reference(r, 6)
        loss = LayeredBandLoss(ref, r.sample_rate, audibility=True)
        rows = []
        args.output.mkdir(parents=True, exist_ok=True)
        trials = [base] + [
            dict(
                base,
                bloom_rate=rate,
                bloom_energy_acceleration=n,
                body_excitation=excitation,
            )
            for rate, n, excitation in product((4, 5.5, 7), (0.05, 0.1, 0.15), (3, 4))
        ]
        if args.excitation_shape:
            trials = [base] + [
                dict(base, body_excitation_centre=centre, bloom_rate=rate)
                for centre, rate in product((700, 1000, 1500), (3, 4, 5.5))
            ]
        for index, p in enumerate(trials):
            db = np.array(
                [loss.envelopes(r.render(p, 6, seed)) for seed in (1675, 1982)]
            )
            offset = ((db - loss.target) * loss.weights).mean(axis=-1)
            # Optimistic high-bar headroom is a coarse feasibility screen only.
            headroom = []
            for low, high in zip(loss.edges[-3:-1], loss.edges[-2:]):
                levels = [
                    p[f"resolved_level_{i}"]
                    for i in range(32)
                    if low <= p[f"resolved_frequency_{i}"] < high
                ]
                headroom.append(6 - min(levels) if levels else 0)
            unavailable = np.maximum(0, -offset[:, -2:] - headroom)
            score = loss.score_db(db, True) + 0.5 * np.sqrt(np.mean(unavailable**2))
            rows.append(
                dict(
                    parameters=p,
                    score=float(score),
                    shape=loss.score_db(db, True),
                    absolute=loss.score_db(db),
                    unavailable_gain_db=unavailable.tolist(),
                )
            )
            print(
                json.dumps(
                    dict(
                        trial=index,
                        score=score,
                        rate=p["bloom_rate"],
                        concentration=p["bloom_energy_acceleration"],
                        excitation=p["body_excitation"],
                    )
                ),
                flush=True,
            )
        rows.sort(key=lambda x: x["score"])
        (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
        results = [
            polish(
                r,
                loss,
                ref,
                row,
                args.output / f"trial-{index}",
                anchors=[120, 360, 900, 2000, 4500, 7000, 10500, 14000],
                fit_seeds=1 if args.excitation_shape else 2,
            )
            for index, row in enumerate(rows[: args.finalists])
        ]
        (args.output / "polished.json").write_text(json.dumps(results, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--finalists", type=int, default=2)
    parser.add_argument("--excitation-shape", action="store_true")
    run(parser.parse_args())
