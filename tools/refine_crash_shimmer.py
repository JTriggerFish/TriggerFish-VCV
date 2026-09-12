"""Saved-crash texture and shared-series refinement using the exact Wasm voice.

No per-mode frequency/decay fitting, EQ optimization, or gain normalization.
Only writes unpublished build artifacts; publication requires a separate review.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from crash_refinement_contracts import training_seeds, audit_seeds, prepare_output

from crash_refinement_search import (
    search_metadata,
    fit_locked_texture,
    fit_shared_series,
    fit_sparse_decay,
    fit_upper_balance,
)
from crash_texture_diagnostics import CrashBalance, plots, spectrograms
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def run(args):
    torch.set_num_threads(1)
    output = args.output
    prepare_output(output)
    fit = json.loads(args.fit.read_text(encoding="utf8"))
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(renderer, fit)
        reference = aligned_reference(renderer, 6)
        objective = CrashBalance(reference, renderer.sample_rate)
        seeds = training_seeds(fit["controls"]["event"]["seed"])
        holdout_seed = audit_seeds(seeds[0])[0]
        rows = []

        def evaluate(p, stage):
            components = [
                objective.components(saved.render(p, 6, seed)) for seed in seeds
            ]
            row = dict(
                stage=stage,
                parameters=p,
                components=components,
                score=float(np.mean([c["score"] for c in components])),
            )
            rows.append(row)
            if len(rows) % 10 == 0:
                (output / "progress.json").write_text(
                    json.dumps(rows, indent=2), encoding="utf8"
                )
            if len(rows) == 1 or row["score"] <= min(r["score"] for r in rows):
                print(
                    json.dumps(
                        dict(
                            stage=stage,
                            evaluation=len(rows),
                            score=row["score"],
                            components=components,
                        )
                    ),
                    flush=True,
                )
            return row["score"]

        evaluate(saved.initial, "baseline")
        baseline = saved.render(saved.initial, 6)
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        write_wav(output / "before.wav", AudioBuffer(baseline, renderer.sample_rate))
        (output / "source.fit.json").write_text(
            json.dumps(fit, indent=2), encoding="utf8"
        )
        if args.low_decay:
            for scale in (0.65, 0.8, 0.9):
                evaluate(
                    dict(
                        saved.initial,
                        body_decay_seconds_0=saved.initial["body_decay_seconds_0"]
                        * scale,
                    ),
                    f"low T60 multiplier {scale}",
                )
        elif args.locked_texture and args.budget:
            fit_locked_texture(saved.initial, evaluate, args.budget)
        elif args.upper_balance:
            fit_upper_balance(saved.initial, evaluate)
        elif args.sparse_decay:
            fit_sparse_decay(saved.initial, rows, evaluate, args.fine_decay)
        elif args.budget:
            fit_shared_series(saved.initial, rows, evaluate, args.budget)
        best = min(rows, key=lambda r: r["score"])
        audio = saved.render(best["parameters"], 6)
        holdout = objective.components(
            saved.render(best["parameters"], 6, holdout_seed)
        )
        report = dict(
            source=fit,
            rows=rows,
            selected=best,
            holdout=holdout,
            holdout_seed=holdout_seed,
            objective=objective.specification,
            training_seeds=seeds,
            **search_metadata(args),
            renderer_sha256=renderer.metadata["rendererSha256"],
        )
        (output / "report.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        (output / "progress.json").write_text(
            json.dumps(rows, indent=2), encoding="utf8"
        )
        (output / "candidate.fit.json").write_text(
            json.dumps(
                saved.snapshot(best["parameters"], args.name),
                indent=2,
            ),
            encoding="utf8",
        )
        write_wav(output / "candidate.wav", AudioBuffer(audio, renderer.sample_rate))
        plots(
            {"Reference": reference, "Before": baseline, "Candidate": audio},
            renderer.sample_rate,
            output,
        )
        spectrograms(reference, audio, renderer.sample_rate, output)
        print(
            json.dumps(
                dict(selected=best["stage"], score=best["score"], holdout=holdout)
            ),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fit", type=Path, default=Path("workbench/web/crash_calibration.fit.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("build/crash-shimmer-refinement")
    )
    parser.add_argument("--budget", type=int, default=180)
    parser.add_argument("--name", default="Crash — clearer rings and shimmer")
    stage = parser.add_mutually_exclusive_group()
    stage.add_argument(
        "--low-decay",
        action="store_true",
        help="Shorten the low endpoint after reviewing an interior knot",
    )
    stage.add_argument(
        "--locked-texture",
        action="store_true",
        help="Only shared pitch/stretch and two T60 endpoints",
    )
    stage.add_argument("--sparse-decay", action="store_true")
    stage.add_argument(
        "--upper-balance",
        action="store_true",
        help="Screen broad upper prominence slopes only",
    )
    stage.add_argument(
        "--fine-decay",
        action="store_true",
        help="Smaller one-knot trials; implies --sparse-decay",
    )
    args = parser.parse_args()
    if args.budget < 0:
        parser.error("--budget must be nonnegative")
    args.sparse_decay |= args.fine_decay
    run(args)
