"""Gentler paired crash: separate beat depth, shared geometry, broad balance.

Uses exact Wasm at the reference gesture. No gain matching, per-mode decay or
independent ridge placement. Proposals are saved, never automatically published.
"""

import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize, lsq_linear

from triggerfish_percussion.coarse_observation_fit import (
    CoarseObservationBasis,
    interpolation_weights,
)
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.low_mode_beating import LowModeBeating
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_fit_baseline import check_reference
from crash_beating_common import CrashObjective
from fit_stretched_gong import checkpoint

KNOTS = (125, 240, 450, 1000, 3500, 15000)


def broad_start(base, middle_scale):
    """One smooth middle-frequency warp; six broad observation coordinates."""
    p = dict(base)
    active = [i for i in range(32) if p[f"resolved_level_{i}"] > -71.99]
    for i in active:
        f = p[f"resolved_frequency_{i}"]
        weight = np.interp(np.log(f), np.log([180, 400, 1200, 4000]), [0, 1, 1, 0])
        p[f"resolved_frequency_{i}"] = float(f * middle_scale**weight)
    weights = interpolation_weights(
        [p[f"resolved_frequency_{i}"] for i in active], KNOTS
    )
    levels = np.array([p[f"resolved_level_{i}"] for i in active])
    coordinates = lsq_linear(weights, 10 ** (levels / 20), bounds=(0.006, 1.99)).x
    for i, level in zip(active, 20 * np.log10(weights @ coordinates)):
        p[f"resolved_level_{i}"] = float(level)
    return p, coordinates


def screen_rings(renderer, objective, beats, base):
    """Screen only pair depth and frequency-rate tilt, with fixed spectrum."""
    rows = []
    # Depth is intentionally bounded below the rejected .5 trial. Beat
    # speed at 125 Hz stays fixed; inspect all four bands, not only the low.
    for depth in (0.1, 0.2, 0.3):
        for tilt in (0, 0.25, 0.5):
            p = dict(base, field_beat_depth=depth, field_beat_rate_tilt=tilt)
            audio = renderer.render(p, 6)
            bands = beats.analyze(audio)
            c = objective.components(audio)
            score = objective.score_components(c) + 0.3 * beats.score_rows(bands)
            rows.append(
                dict(
                    parameters=p,
                    score=score,
                    components=c,
                    depth=depth,
                    tilt=tilt,
                    beating=beats.score_rows(bands),
                )
            )
            print(
                json.dumps({k: v for k, v in rows[-1].items() if k != "parameters"}),
                flush=True,
            )
    return rows


def balance_trial(renderer, objective, beats, base, seeds, scale, output, rows):
    """Optimize six observation weights for one smoothly warped geometry."""
    p, coordinates = broad_start(base, scale)
    basis = CoarseObservationBasis(renderer, p, 6, seeds, KNOTS)
    tensors = [(a, c) for a, c in basis.bases.values()]
    evaluations = []

    def evaluate(db):
        amplitudes = 10 ** (np.array(db) / 20)
        audio = [a + (amplitudes - basis.amplitudes) @ c for a, c in tensors]
        # Full and attack spectra drive balance, with a small envelope
        # term. Texture is audited after fitting, not chased by levels.
        parts = [
            dict(
                mel=objective.mel.score(a),
                attack=objective.attack.score(a[: objective.attack_frames]),
                bloom=float(np.linalg.norm(objective.shape.residual(a))),
            )
            for a in audio
        ]
        value = float(
            np.mean([r["mel"] + 0.3 * r["attack"] + 0.05 * r["bloom"] for r in parts])
        )
        evaluations.append(dict(value=value, db=list(db), components=parts))
        return value

    minimize(
        evaluate,
        20 * np.log10(coordinates),
        method="Powell",
        bounds=[(-45, 5.99)] * len(KNOTS),
        options=dict(maxfev=160, xtol=0.025, ftol=0.001),
    )
    chosen = min(evaluations, key=lambda r: r["value"])
    fitted = basis.parameters(10 ** (np.array(chosen["db"]) / 20))
    audio = renderer.render(fitted, 6)
    predicted = (
        tensors[0][0]
        + (10 ** (np.array(chosen["db"]) / 20) - basis.amplitudes) @ tensors[0][1]
    )
    error = np.linalg.norm(audio - predicted) / max(np.linalg.norm(audio), 1e-12)
    if error > 3e-4:
        raise ValueError(f"Observation basis render mismatch: {error}")
    directory = output / f"middle-{scale}"
    checkpoint(
        renderer,
        objective,
        directory,
        "Crash — gentle ring and fuller body",
        fitted,
        reference,
        [
            dict(stage="depth/tilt screen", trials=rows),
            dict(
                stage="broad balance",
                middle_scale=scale,
                knots=KNOTS,
                basis_validation=basis.validation,
                actual_error=float(error),
                evaluations=evaluations,
                selected=chosen,
            ),
        ],
    )
    summary = dict(
        directory=str(directory),
        scale=scale,
        score=chosen["value"],
        components=objective.components(audio),
        beating=beats.score(audio),
    )
    print(json.dumps(summary), flush=True)
    return summary


def run():
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    output = Path("build/crash-gentle-ring")
    try:
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(reference, renderer.sample_rate, 60)
        beats = LowModeBeating(reference, renderer.sample_rate)
        source_path = output / "before" / "search.json"
        if source_path.exists():
            source = json.loads(source_path.read_text(encoding="utf8"))
            check_reference(source["metadata"], renderer.metadata)
            base = source["parameters"]
            if set(base) != set(renderer.initial):
                raise ValueError(
                    "Archived starting point has a different parameter surface"
                )
        else:
            base = dict(renderer.initial)
        checkpoint(
            renderer,
            objective,
            output / "before",
            "Previous paired crash",
            base,
            reference,
            [],
        )
        rows = screen_rings(renderer, objective, beats, base)
        best = min(rows, key=lambda r: r["score"])
        seed = renderer.metadata["event"]["seed"]
        seeds = (seed, (seed + 307) & 0xFFFFFFFF)
        trials = []
        for scale in (0.9, 1, 1.1):
            trials.append(
                balance_trial(
                    renderer,
                    objective,
                    beats,
                    best["parameters"],
                    seeds,
                    scale,
                    output,
                    rows,
                )
            )
        (output / "summary.json").write_text(json.dumps(trials, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    run()
