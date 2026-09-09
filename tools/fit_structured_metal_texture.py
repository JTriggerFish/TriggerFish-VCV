"""Fit whole modal families and global texture, never individual upper ridges.

Every comparison renders the current Wasm at the fixed reference gesture/gain.
Artifacts are candidates only; publishing requires independent review.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


class Objective:
    units = "bloom/10 + 0.6 texture + 0.15 Mel"

    def __init__(self, reference, rate):
        self.shape = SpectralBloomLoss(reference, rate)
        self.texture = ModalTextureLoss(reference, rate)
        self.mel = AuralossMel(reference, rate)
        self.specification = dict(
            version="structured-metal-texture-v1",
            bloom=self.shape.specification,
            texture=self.texture.specification,
            mel=self.mel.specification,
            weights=[0.1, 0.6, 0.15],
        )

    def components(self, audio):
        return dict(
            bloom=float(np.linalg.norm(self.shape.residual(audio))),
            texture=self.texture.score(audio),
            mel=self.mel.score(audio),
        )

    def score(self, audio):
        c = self.components(audio)
        return 0.1 * c["bloom"] + 0.6 * c["texture"] + 0.15 * c["mel"]

    def residual(self, audio, regions=None):
        return np.array([np.sqrt(self.score(audio))])


def generated(renderer, base, settings):
    points = renderer.request(
        command="modalTemplate",
        settings=dict(
            settings, family="harmonic", minimumFrequency=1, maximumFrequency=15000
        ),
    )["points"]
    result = dict(base)
    for i in range(32):
        result[f"resolved_level_{i}"] = points[i]["level"] if i < len(points) else -72
        result[f"resolved_turbulence_{i}"] = 1
        result[f"resolved_allocation_{i}"] = 1
        if i < len(points):
            result[f"resolved_frequency_{i}"] = points[i]["frequency"]
    return result


def screen_series(renderer, base, objective, save):
    """Compare whole families; the old irregular fit remains a separate baseline."""
    candidates = []
    for root in (120, 128, 137):
        for count, stretch in ((24, 0.7), (32, 0.4)):
            settings = dict(
                fundamental=root,
                count=count,
                stretch=stretch,
                harmonicCore=4,
                level=-12,
                rolloff=0,
            )
            p = generated(renderer, base, settings)
            p.update(
                field_phase_bandwidth=0.02,
                field_turbulence=0.8,
                field_distribution=2,
                field_packet_spread=2,
                field_satellite_density=0.5,
            )
            trial = save(
                f"series-{root}-{count}",
                p,
                [dict(stage="whole series", settings=settings)],
            )
            trial.loss = objective.shape
            polish_coarse(trial, (120, 600, 3000, 15000))
            trial.loss = objective
            value = objective.score(trial.audio(trial.parameters))
            print(
                json.dumps(dict(stage="series", settings=settings, score=value)),
                flush=True,
            )
            candidates.append((value, trial))
    return min(candidates, key=lambda row: row[0])[1]


def screen_texture(renderer, initial, objective, save):
    """Change texture only; never overwrite the starting trial's parameter object."""
    best = initial
    best_score = objective.score(best.audio(best.parameters))
    base, history = dict(initial.parameters), list(initial.history)
    for layout in (0, 1, 2):
        for density in (0.15, 0.45, 0.85):
            for blur in (0, 0.01, 0.04, 0.16):
                p = dict(
                    base,
                    field_distribution=layout,
                    field_satellite_density=density,
                    field_phase_bandwidth=blur,
                )
                value = objective.score(renderer.render(p, 6))
                record = dict(
                    stage="texture screen",
                    layout=layout,
                    density=density,
                    blur=blur,
                    score=value,
                )
                print(json.dumps(record), flush=True)
                if value < best_score:
                    best = save("texture-best", p, history + [record])
                    best_score = value
    return best


def refine_shared(search, objective, budget):
    """Optimize existing global controls, followed by a guarded broad level fit."""
    bounds = dict(
        bloom_rate=(0.1, 16),
        bloom_energy_acceleration=(0, 0.25),
        body_brightness=(-60, 6),
        body_excitation_centre=(150, 4000),
        field_turbulence=(0.1, 2),
        field_turbulence_slope=(-0.2, 0.8),
        field_packet_spread=(0.1, 6),
        field_phase_bandwidth=(0, 0.3),
        body_decay_seconds_0=(0.3, 30),
        body_decay_seconds_7=(0.1, 15),
    )
    bounds = {
        k: (min(lo, search.parameters[k]), max(hi, search.parameters[k]))
        for k, (lo, hi) in bounds.items()
    }
    refine_scalar(search, bounds, budget=budget, step=0.004, method="Powell")
    before = dict(search.parameters)
    before_score = objective.score(search.audio(before))
    search.loss = objective.shape
    polish_coarse(search, (120, 600, 3000, 15000))
    search.loss = objective
    proposal_score = objective.score(search.audio(search.parameters))
    if proposal_score > before_score:
        search.parameters = before
    search.history.append(
        dict(
            stage="broad prominence acceptance",
            before=before_score,
            proposal=proposal_score,
            accepted=proposal_score <= before_score,
        )
    )
    search.save()


def compare_seeds(renderer, objective, base, candidate):
    """Report components separately; a score is not listening approval."""
    rows = []
    for offset in (0, 307, 911):
        seed = (renderer.metadata["event"]["seed"] + offset) & 0xFFFFFFFF
        rows.append(
            dict(
                seed=seed,
                baseline=objective.components(renderer.render(base, 6, seed)),
                candidate=objective.components(renderer.render(candidate, 6, seed)),
            )
        )
    return dict(
        objective=objective.specification, comparison=rows, listening_approved=False
    )


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        reference = aligned_reference(renderer, 6)
        objective = Objective(reference, renderer.sample_rate)
        base = dict(renderer.initial)
        if args.fit:
            source = json.loads(args.fit.read_text(encoding="utf8"))
            if (
                source["reference"]["sha256"]
                != renderer.metadata["reference"]["sha256"]
            ):
                raise ValueError("User fit and reference target differ")
            base.update(
                {
                    k: v
                    for n in source["instrument"]["nodes"]
                    for k, v in n["parameters"].items()
                }
            )

        def save(label, params, history):
            return checkpoint(
                renderer,
                objective,
                args.output / label,
                args.target.title() + " — structured beating trial",
                params,
                reference,
                history,
            )

        baseline = save(
            "baseline",
            base,
            [dict(stage="current-engine starting point", source=str(args.fit))],
        )
        print(
            json.dumps(
                dict(
                    stage="baseline",
                    components=objective.components(baseline.audio(base)),
                )
            ),
            flush=True,
        )
        initial = (
            screen_series(renderer, base, objective, save)
            if args.target == "crash"
            else baseline
        )
        best = screen_texture(renderer, initial, objective, save)
        search = save("candidate", best.parameters, best.history)
        refine_shared(search, objective, args.budget)
        verify_candidate(renderer, search.output)
        report = compare_seeds(renderer, objective, base, search.parameters)
        (args.output / "comparison.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(json.dumps(report), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("--fit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=260)
    run(parser.parse_args())
