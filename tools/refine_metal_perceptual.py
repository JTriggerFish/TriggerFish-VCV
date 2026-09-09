"""Refine an exact-render metallic checkpoint with the library mel objective.

No source gains, event controls, frequencies, knot counts or per-mode damping
change. Positive observation amplitudes use analysis autograd; texture/contact
use bounded Powell searches over the actual Wasm, not a surrogate instrument.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def scalar_stage(search, name, bounds, budget):
    """Normalized coordinates, exact objective, best-evaluated checkpoint."""
    descriptors = {d["key"]: d for d in search.renderer.metadata["descriptors"]}
    for key, (low, high) in bounds.items():
        d = descriptors[key]
        if not d["minimum"] <= low < high <= d["maximum"]:
            raise ValueError(f"Search outside exposed control range: {key}")
    keys = list(bounds)
    low, high = np.array(list(bounds.values())).T
    seed = dict(search.parameters)
    initial = np.array([seed[k] for k in keys])
    start_score = search.loss.score(search.audio(seed))
    if not np.isfinite(start_score):
        raise ValueError("Nonfinite scalar baseline")
    best_score, best_parameters = start_score, seed
    evaluations = 0

    def objective(unit):
        nonlocal best_score, best_parameters, evaluations
        parameters = dict(seed, **dict(zip(keys, (low + unit * (high - low)).tolist())))
        score = search.loss.score(search.audio(parameters))
        if not np.isfinite(score):
            raise ValueError("Nonfinite scalar candidate")
        evaluations += 1
        if score < best_score:
            best_score, best_parameters = score, parameters
        if evaluations % 10 == 0:
            progress = dict(
                stage=name,
                evaluations=evaluations,
                before=start_score,
                best=best_score,
                parameters=best_parameters,
            )
            (search.output / "scalar-progress.json").write_text(json.dumps(progress))
            print(
                json.dumps({k: v for k, v in progress.items() if k != "parameters"}),
                flush=True,
            )
        return score

    start = np.clip((initial - low) / (high - low), 0, 1)
    influence = []
    for index, key in enumerate(keys):
        minus, plus = start.copy(), start.copy()
        minus[index] = max(0, start[index] - 0.02)
        plus[index] = min(1, start[index] + 0.02)
        scores = [objective(point) for point in (minus, plus)]
        influence.append(
            dict(
                parameter=key,
                step=(plus[index] - minus[index]) * (high[index] - low[index]) / 2,
                minus=scores[0],
                plus=scores[1],
                baseline=start_score,
            )
        )
    result = minimize(
        objective,
        start,
        method="Powell",
        bounds=[(0, 1)] * len(keys),
        options=dict(maxfev=budget, xtol=0.005, ftol=0.0005),
    )
    search.parameters = best_parameters
    search.history.append(
        dict(
            stage=name,
            method="bounded Powell / exact renderer",
            bounds=bounds,
            fixed_parameters={k: v for k, v in seed.items() if k not in bounds},
            before=start_score,
            after=best_score,
            evaluations=evaluations,
            influence=influence,
            solver_message=str(result.message),
            objective=search.loss.specification,
        )
    )
    search.save()


def checkpoint(search, name, baseline, reference):
    root = search.output
    search.output = root / name
    search.output.mkdir(parents=True, exist_ok=True)
    search.save()
    write_wav(
        search.output / "reference.wav",
        AudioBuffer(reference, search.renderer.sample_rate),
    )
    verify_candidate(search.renderer, search.output)
    balance = MetallicBalanceLoss(reference, search.renderer.sample_rate, "erb", True)
    decay = BandDecayShapeLoss(reference, search.renderer.sample_rate)
    report = {}
    for label, parameters in [("baseline", baseline), ("candidate", search.parameters)]:
        audio = search.audio(parameters)
        report[label] = dict(
            mel=search.loss.score(audio),
            balance=balance.diagnostics(audio),
            decay=decay.diagnostics(audio),
        )
    (search.output / "review.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(dict(checkpoint=name, **report)), flush=True)
    search.output = root


def fit(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.resume)
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        search = Search(
            renderer,
            AuralossMel(reference, renderer.sample_rate),
            args.output,
            seconds,
            f"{args.target.title()} — relaxed reference refinement",
            (None,),
        )
        args.output.mkdir(parents=True, exist_ok=True)
        search.parameters = dict(saved["parameters"])
        if args.cascade_min:
            search.parameters["bloom_rate"] = max(
                args.cascade_min, search.parameters["bloom_rate"]
            )
        baseline = dict(renderer.initial)
        search.history.append(dict(parent=str(args.resume.resolve())))
        polish_observation_autograd(search, iterations=100)
        checkpoint(search, "mel-observation", baseline, reference)
        for turn in range(args.rounds):
            # Centre stays fixed: otherwise it is redundant with level.
            scalar_stage(
                search,
                f"round {turn+1}: texture",
                dict(
                    field_turbulence=(0.05, 4),
                    field_turbulence_slope=(0, 1),
                    field_phase_bandwidth=(0, 1.2),
                    field_packet_spread=(0.2, 8),
                ),
                160,
            )
            scalar_stage(
                search,
                f"round {turn+1}: transport",
                dict(
                    bloom_rate=(max(0.01, args.cascade_min), 8),
                    body_brightness=(-36, 18),
                ),
                100,
            )
            decay = {
                "body_decay_seconds_0": (0.2, 30),
                "body_decay_seconds_7": (0.15, 8),
            }
            for knot in range(1, 7):
                if search.parameters[f"body_decay_active_{knot}"] >= 0.5:
                    decay[f"body_decay_seconds_{knot}"] = (0.2, 25)
            scalar_stage(search, f"round {turn+1}: damping", decay, 100)
            scalar_stage(
                search,
                f"round {turn+1}: contact",
                dict(
                    impact_tone_noise=(0, 1),
                    impact_width=(0.25, 4),
                    impact_noise_tilt=(-18, 18),
                    impact_chirp_pitch=(0.1, 4),
                ),
                100,
            )
            polish_observation_autograd(search, iterations=100)
            checkpoint(search, f"round-{turn+1}", baseline, reference)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--cascade-min", type=float, default=0)
    args = parser.parse_args()
    if not 0 <= args.cascade_min < 8:
        parser.error("Cascade minimum must be in [0, 8)")
    fit(args)
