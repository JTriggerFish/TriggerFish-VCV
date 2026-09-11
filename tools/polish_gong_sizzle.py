"""Retain metallic spectral detail while fitting upper colour and modulation.

One smooth observation adjustment, not independent modal placement. Preserve
the accepted low body and reject loss of > .75 dB fine spectral contrast in
either upper band relative to the incumbent. These are explicit design guards,
not perceptual thresholds. The fitting objective remains UpperSizzleLoss.
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
from audit_gong_sizzle import ridge_contrast
from fit_stretched_gong import checkpoint


def colour(base, middle_db):
    """Smooth observation curve: unchanged below 900 Hz and above 12 kHz."""
    result = dict(base)
    knots = np.log([900, 4000, 12000])
    for i in range(32):
        key = f"resolved_level_{i}"
        f = np.log(base[f"resolved_frequency_{i}"])
        j = np.clip(np.searchsorted(knots, f) - 1, 0, 1)
        x = np.clip((f - knots[j]) / (knots[j + 1] - knots[j]), 0, 1)
        x = x * x * (3 - 2 * x)
        gain = middle_db * (x if j == 0 else 1 - x)
        if base[key] > -71.99:
            result[key] = float(np.clip(base[key] + gain, -72, 6))
    return result


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = json.loads(args.baseline.read_text())
        ref = aligned_reference(r, 6)
        loss = UpperSizzleLoss(ref, r.sample_rate)
        seeds = [1675, 1982, 2586]
        original_audio = [r.render(base, 6, s) for s in seeds]
        originals = [loss.measure(a) for a in original_audio]
        contrast = np.array([ridge_contrast(a, r.sample_rate) for a in original_audio])
        args.output.mkdir(parents=True, exist_ok=True)
        rows = []

        def evaluate(name, p):
            audio = [r.render(p, 6, s) for s in seeds]
            metrics = [loss.diagnostics(a, b) for a, b in zip(audio, originals)]
            ridge = np.array([ridge_contrast(a, r.sample_rate) for a in audio])
            scores = [m["score"] for m in metrics]
            eligible = max(m["low_max_change_db"] for m in metrics) <= 2
            eligible &= np.min(ridge - contrast) >= -0.75
            row = dict(
                name=name,
                parameters=p,
                metrics=metrics,
                ridge_contrast_db=ridge.tolist(),
                eligible=bool(eligible),
                score=float(np.mean(scores) + 0.5 * max(scores)),
            )
            rows.append(row)
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            print(
                json.dumps(
                    dict(name=name, score=row["score"], eligible=row["eligible"])
                ),
                flush=True,
            )

        evaluate("baseline", base)
        for blur, tilt, middle in product(
            [0.01, 0.035, 0.1], [0.5, 1, 1.5, 2], [-6, -3, 0]
        ):
            p = colour(base, middle)
            p.update(field_phase_bandwidth=blur, field_phase_tilt=tilt)
            evaluate(f"blur {blur}, tilt {tilt}, middle {middle} dB", p)
        best = min(
            (row for row in rows if row["eligible"]), key=lambda row: row["score"]
        )
        checkpoint(
            r,
            loss,
            args.output / "candidate",
            "Gong — brighter sizzle",
            best["parameters"],
            ref,
            [
                dict(
                    stage="guarded upper colour and texture",
                    source=str(args.baseline),
                    objective=loss.specification,
                    seeds=seeds,
                    ridge_loss_guard_db=0.75,
                    curve_hz=[900, 4000, 12000],
                    selected=best["name"],
                    baseline_score=rows[0]["score"],
                    score=best["score"],
                )
            ],
        )
        verify_candidate(r, args.output / "candidate")
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
