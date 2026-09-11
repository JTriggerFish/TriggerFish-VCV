"""Compare bounded ridge movement and blur on the exact gong, without publishing.

Frequencies, levels, damping, excitation and bloom remain fixed. This probe is
a controlled model experiment, not a claim of a completed perceptual fit.
"""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.upper_sizzle_loss import UpperSizzleLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from audit_gong_sizzle import brightness, ridge_contrast
from fit_stretched_gong import checkpoint


def candidates(base, blurs):
    yield "accepted", dict(base)
    yield "no blur or movement", dict(base, field_phase_bandwidth=0)
    for blur, depth, rate, sharing in product(
        blurs, [0.5, 1, 1.5], [10, 40, 100], [0.25, 0.75]
    ):
        yield f"blur {blur}, movement {depth}, speed {rate}, sharing {sharing}", dict(
            base,
            field_phase_bandwidth=blur,
            field_motion_depth=depth,
            field_motion_rate=rate,
            field_motion_sharing=sharing,
        )


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        base = dict(r.initial)
        ref = aligned_reference(r, 6)
        loss = UpperSizzleLoss(ref, r.sample_rate)
        seeds = [1675, 1982, 2586]
        originals = [loss.measure(r.render(base, 6, s)) for s in seeds]
        checkpoint(
            r, loss, args.output / "baseline", "Gong — accepted baseline", base, ref, []
        )
        rows = []
        for name, p in candidates(base, args.blurs):
            measurements = []
            for s, original in zip(seeds, originals):
                response = r.request(parameters=p, seconds=6, seed=s)
                audio = r.decode(response["pcm"])
                m = loss.diagnostics(audio, original)
                m.update(
                    brightness=brightness(audio, r.sample_rate)[0],
                    ridge_contrast_db=ridge_contrast(audio, r.sample_rate),
                    render_ms=response["elapsedMs"],
                    peak_db=float(20 * np.log10(max(abs(audio)))),
                )
                measurements.append(m)
            scores = [m["score"] for m in measurements]
            row = dict(
                name=name,
                parameters=p,
                metrics=measurements,
                score=float(np.mean(scores) + 0.5 * max(scores)),
                eligible=max(m["low_max_change_db"] for m in measurements) <= 2,
            )
            # Protect the existing metallic fine detail, not only broad colour.
            if rows:
                original_ridges = np.array(
                    [m["ridge_contrast_db"] for m in rows[0]["metrics"]]
                )
                row["eligible"] &= bool(
                    np.min(
                        np.array([m["ridge_contrast_db"] for m in measurements])
                        - original_ridges
                    )
                    >= -0.75
                )
            rows.append(row)
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            print(
                json.dumps(
                    dict(name=name, score=row["score"], eligible=row["eligible"])
                ),
                flush=True,
            )
        selected = min(
            (row for row in rows if row["eligible"]), key=lambda row: row["score"]
        )
        checkpoint(
            r,
            loss,
            args.output / "candidate",
            "Gong — ridge movement experiment",
            selected["parameters"],
            ref,
            [
                dict(
                    stage="bounded motion probe",
                    seeds=seeds,
                    selected=selected["name"],
                    objective=loss.specification,
                    ridge_guard_db=0.75,
                    score=selected["score"],
                )
            ],
        )
        verify_candidate(r, args.output / "candidate")
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blurs", type=float, nargs="+", default=[0, 0.01, 0.035])
    run(parser.parse_args())
