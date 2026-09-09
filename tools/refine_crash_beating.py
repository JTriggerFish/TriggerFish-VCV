"""Crash-only structured search: low blur, whole series, fixed levels and gesture.

The gong and DSP are untouched. Search stages are independently resumable;
checkpoints are proposals, never automatic workbench publication.
"""

import argparse
import json
import os
from pathlib import Path

import torch

from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.coarse_observation_fit import (
    polish_coarse,
)
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from fit_structured_metal_texture import generated
from refine_coarse_metal import damping_variant

from crash_beating_common import CrashObjective, KNOTS, project_levels


def screen_texture(renderer, base, objective, save):
    rows = []
    for layout in (0, 1, 2):
        for density in (0.45, 0.85, 1):
            for blur in (0, 0.001, 0.004, 0.01):
                p = dict(
                    base,
                    field_distribution=layout,
                    field_satellite_density=density,
                    field_phase_bandwidth=blur,
                )
                audio = renderer.render(p, 6)
                row = dict(
                    layout=layout,
                    density=density,
                    blur=blur,
                    score=objective.score(audio),
                    components=objective.components(audio),
                    parameters=p,
                )
                rows.append(row)
                print(
                    json.dumps({k: v for k, v in row.items() if k != "parameters"}),
                    flush=True,
                )
    best = min(rows, key=lambda row: row["score"])
    save(
        "texture",
        best["parameters"],
        [dict(stage="fixed-geometry low-blur screen", trials=rows)],
    )
    return best["parameters"]


def screen_families(renderer, base, objective, save):
    rows = []
    for root in (110, 120, 130, 145):
        for core in (2, 4, 6):
            settings = dict(fundamental=root, count=32, harmonicCore=core)
            settings["stretch"] = renderer.request(
                command="modalTemplateStretch",
                settings=dict(settings, topFrequency=13500),
            )["stretch"]
            p = project_levels(
                generated(renderer, base, dict(settings, level=-12, rolloff=0)), base
            )
            trial = save(
                f"series-{root}-{core}",
                p,
                [dict(stage="whole series", settings=settings)],
            )
            trial.loss = objective.shape
            polish_coarse(trial, KNOTS)
            trial.loss = objective
            audio = renderer.render(trial.parameters, 6)
            row = dict(
                settings=settings,
                score=objective.score(audio),
                components=objective.components(audio),
                directory=str(trial.output),
            )
            rows.append(row)
            print(json.dumps(row), flush=True)
    return sorted(rows, key=lambda row: row["score"])


def refine(search, objective, budget, stable=False):
    """Low blur is a hypothesis constraint; all other bounds are explicit."""
    search.loss = objective
    bounds = dict(
        bloom_rate=(0.1, 16),
        bloom_energy_acceleration=(0, 0.25),
        body_brightness=(-48, 12),
        body_excitation_centre=(150, 4500),
        field_turbulence=(0.1, 2.5),
        field_turbulence_slope=(0, 1),
        field_packet_spread=(0.1, 5),
        field_phase_bandwidth=(0, 0.01),
        body_decay_seconds_0=(0.3, 30),
        body_decay_seconds_7=(0.1, 10),
    )
    if stable:
        bounds.pop("field_phase_bandwidth")
        bounds["bloom_energy_acceleration"] = (0, 1)
    # Preserve a valid starting point; names/ranges are checked by ParameterBox.
    refine_scalar(search, bounds, budget=budget, step=0.004, method="Powell")
    guarded_polish(search, objective)


def guarded_polish(search, objective):
    """Keep a smooth observation edit only when the full proposal score improves."""
    before = dict(search.parameters)
    initial = objective.score(search.audio(before))
    search.loss = objective.shape
    polish_coarse_mel(search, objective.mel, KNOTS)
    proposal = objective.score(search.audio(search.parameters))
    if proposal > initial:
        search.parameters = before
    search.loss = objective
    search.history.append(
        dict(
            stage="guarded spectral observation polish",
            before=initial,
            proposal=proposal,
            accepted=proposal <= initial,
        )
    )
    search.save()


def refine_contact(search, objective, budget):
    """Revisit the shared exciter after body fitting, without changing gesture."""
    search.loss = objective
    refine_scalar(
        search,
        dict(
            impact_tone_noise=(0, 1),
            impact_width=(0.25, 2),
            impact_noise_tilt=(-18, 18),
            impact_chirp_pitch=(0.25, 4),
            direct_gain=(0, 0.5),
        ),
        budget=budget,
        step=0.004,
        method="Powell",
    )
    guarded_polish(search, objective)


def refine_decay(search, objective, budget):
    """One shared interior knot, only after the two-endpoint body/contact fit."""
    original = dict(search.parameters)
    before = objective.score(search.audio(original))
    search.parameters = damping_variant(original, [600])
    search.loss = objective
    refine_scalar(
        search,
        dict(
            body_decay_seconds_0=(0.3, 30),
            body_decay_seconds_1=(0.1, 30),
            body_decay_seconds_7=(0.1, 15),
        ),
        budget=budget,
        step=0.004,
        method="Powell",
    )
    guarded_polish(search, objective)
    after = objective.score(search.audio(search.parameters))
    accepted = after < before - 0.02
    if not accepted:
        search.parameters = original
    search.history.append(
        dict(
            stage="extra shared decay knot acceptance",
            frequency=600,
            before=before,
            after=after,
            accepted=accepted,
            minimum_improvement=0.02,
        )
    )
    search.save()


def refine_dynamics(search, objective, budget):
    """Fit development/decay while holding the selected packet texture fixed."""
    search.loss = objective
    refine_scalar(
        search,
        dict(
            bloom_rate=(0.1, 16),
            bloom_energy_acceleration=(0, 1),
            body_brightness=(-48, 12),
            body_excitation_centre=(150, 4500),
            body_decay_seconds_0=(0.3, 30),
            body_decay_seconds_7=(0.1, 10),
        ),
        budget=budget,
        step=0.004,
        method="Powell",
    )
    guarded_polish(search, objective)


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(
            reference, renderer.sample_rate, args.reference_floor_db
        )

        def save(label, p, history):
            return checkpoint(
                renderer,
                objective,
                args.output / label,
                "Crash — low-blur structured trial",
                p,
                reference,
                history,
            )

        if args.stage == "screen":
            base = dict(renderer.initial)
            save("baseline", base, [dict(stage="published starting point")])
            texture = screen_texture(renderer, base, objective, save)
            rows = screen_families(renderer, texture, objective, save)
            (args.output / "series-screen.json").write_text(json.dumps(rows, indent=2))
        else:
            if not args.source:
                raise ValueError("--source checkpoint is required for refinement")
            source = verify_candidate(renderer, args.source)
            search = save("candidate", source["parameters"], source["history"])
            if args.phase_blur is not None:
                if not 0 <= args.phase_blur <= 4:
                    raise ValueError("Phase blur must be in [0, 4]")
                search.parameters["field_phase_bandwidth"] = args.phase_blur
                search.history.append(
                    dict(
                        stage="explicit fixed blur starting point",
                        field_phase_bandwidth=args.phase_blur,
                    )
                )
            if args.stage == "stable":
                search.parameters.update(
                    field_phase_bandwidth=0,
                    field_satellite_density=1,
                    field_distribution=2,
                )
                search.history.append(
                    dict(
                        stage="stable beating hypothesis",
                        fixed=dict(
                            field_phase_bandwidth=0,
                            field_satellite_density=1,
                            field_distribution=2,
                        ),
                    )
                )
                search.loss = objective.shape
                polish_coarse(search, KNOTS)
                refine(search, objective, args.budget, stable=True)
            else:
                action = {
                    "contact": refine_contact,
                    "decay": refine_decay,
                    "dynamics": refine_dynamics,
                    "polish": lambda search, objective, budget: guarded_polish(
                        search, objective
                    ),
                }.get(args.stage, refine)
                action(search, objective, args.budget)
            verify_candidate(renderer, search.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=[
            "screen",
            "refine",
            "contact",
            "stable",
            "decay",
            "dynamics",
            "polish",
        ],
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--budget", type=int, default=420)
    parser.add_argument(
        "--phase-blur",
        type=float,
        help="Explicit experiment starting value; contact/decay stages keep it fixed",
    )
    parser.add_argument(
        "--reference-floor-db",
        type=float,
        help="Reference-derived spectral comparison range; does not alter audio",
    )
    run(parser.parse_args())
