"""Staged exact-Wasm crash/gong fitting; never publishes factory presets.

Shared objective and algorithm, per-target starting patch. Structural texture
is screened before finite-difference refinement; observation uses a validated
linear render basis. Every checkpoint retains the full patch and fixed event.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.stats import qmc
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def checkpoint(search, name, baseline, mel):
    previous = search.output
    search.output = previous / name
    search.output.mkdir(parents=True, exist_ok=True)
    search.save()
    write_wav(
        search.output / "reference.wav",
        AudioBuffer(
            aligned_reference(search.renderer, search.seconds),
            search.renderer.sample_rate,
        ),
    )
    verify_candidate(search.renderer, search.output)
    report = {"stage": name, "evaluations": search.evaluations}
    for label, parameters in [("baseline", baseline), ("candidate", search.parameters)]:
        audio = search.audio(parameters)
        report[label] = dict(
            balance=search.loss.diagnostics(audio), mel=mel.score(audio)
        )
    (search.output / "review.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)
    search.output = previous


def fit(target, output, rounds, resume, cascade_min=0):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{target}-standard", Path.cwd()
    )
    try:
        seconds = 10 if target == "crash" else 8
        reference = aligned_reference(renderer, seconds)
        loss = MetallicBalanceLoss(reference, renderer.sample_rate, "erb", True)
        mel = AuralossMel(reference, renderer.sample_rate)
        output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            loss,
            output,
            seconds,
            f"{target.title()} — relaxed turbulence refit",
            (None,),
        )
        baseline = dict(renderer.initial)
        if resume:
            search.parameters = dict(verify_candidate(renderer, resume)["parameters"])
        if cascade_min:
            search.parameters["bloom_rate"] = max(
                cascade_min, search.parameters["bloom_rate"]
            )
        checkpoint(search, "start", baseline, mel)
        texture = dict(
            field_turbulence=(0.15, 2.5),
            field_turbulence_slope=(0, 1),
            field_phase_bandwidth=(0, 1.2),
            field_packet_spread=(0.5, 8),
        )
        # Noisiness is a level at 1 kHz plus slope; no redundant centre fit.
        if cascade_min:
            # A user-auditioned character constraint, not inferred from the loss.
            # Compensate the excitation distribution before local refinement.
            seed = dict(search.parameters)
            candidates = [
                (
                    f"cascade-{rate}-tilt-{tilt}-accel-{accel}",
                    dict(
                        seed,
                        bloom_rate=rate,
                        body_brightness=tilt,
                        bloom_energy_acceleration=accel,
                    ),
                )
                for rate in (
                    cascade_min,
                    min(8, cascade_min * 1.4),
                    min(8, cascade_min * 2),
                )
                for tilt in (-24, -12, -6, 0)
                for accel in (seed["bloom_energy_acceleration"],)
            ]
            search.screen_candidates(
                "audition-constrained fast cascade starts", candidates
            )
        if not resume:
            seed = dict(search.parameters)
            candidates = []
            for index, unit in enumerate(
                qmc.LatinHypercube(len(texture), seed=731).random(32)
            ):
                values = dict(seed)
                for (key, (low, high)), value in zip(texture.items(), unit):
                    values[key] = float(low + (high - low) * value)
                candidates.append((f"texture-screen-{index}", values))
            contrast = (0.7, 0.65, 2500) if target == "crash" else (0.4, 1, 2000)
            candidates.append(
                (
                    "audition-contrast",
                    dict(
                        seed,
                        **dict(
                            zip(
                                (
                                    "field_turbulence",
                                    "field_turbulence_slope",
                                ),
                                contrast,
                            )
                        ),
                    ),
                )
            )
            search.screen_candidates("relaxed texture starts", candidates)
        for iteration in range(rounds):
            # Explicit broad steps are recorded by Search; low-influence knobs
            # are frozen, not allowed to wander at a numerical boundary.
            search.stage(
                f"pass {iteration+1}: texture",
                texture,
                14,
                difference_step=0.015,
                influence_threshold=0.01,
            )
            motion = dict(
                bloom_rate=(max(0.02, cascade_min), 8),
                body_brightness=(-36, 18),
                body_excitation_centre=(100, 5000),
            )
            search.stage(
                f"pass {iteration+1}: excitation and transport",
                motion,
                18,
                difference_step=0.008,
                influence_threshold=0.01,
            )
            polish_observation_autograd(search, iterations=50)
            decay = {
                "body_decay_seconds_0": (0.2, 30),
                "body_decay_seconds_7": (0.15, 8),
            }
            for knot in range(1, 7):
                if search.parameters[f"body_decay_active_{knot}"] >= 0.5:
                    decay[f"body_decay_seconds_{knot}"] = (0.2, 25)
            search.stage(
                f"pass {iteration+1}: shared damping",
                decay,
                18,
                difference_step=0.005,
                influence_threshold=0.01,
            )
            polish_observation_autograd(search, iterations=40)
            checkpoint(search, f"pass-{iteration+1}", baseline, mel)
        search.save()
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--cascade-min", type=float, default=0)
    args = parser.parse_args()
    if not 0 <= args.cascade_min < 8:
        parser.error("Cascade minimum must be in [0, 8)")
    fit(args.target, args.output, args.rounds, args.resume, args.cascade_min)
