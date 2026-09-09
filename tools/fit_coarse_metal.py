"""Whole stretched-series layouts and shared dynamics, never individual modes.

Archive the current preset; search exact Wasm renders with fixed reference,
gesture and gains. Fit a few broad observation amplitudes, then shared dynamics
and two damping endpoints. Final choice remains a candidate until reviewed.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.coarse_observation_fit import polish_coarse
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel as perceptual_polish


def layouts(
    renderer,
    roots,
    neutral=False,
    original=None,
    counts=(24, 32),
    stretches=(0.25, 0.5, 0.75),
):
    for root in roots:
        for count in counts:
            for stretch in stretches:
                settings = dict(
                    family="harmonic",
                    fundamental=root,
                    count=count,
                    stretch=stretch,
                    harmonicCore=4,
                    level=-20,
                    rolloff=0,
                )
                try:
                    points = renderer.request(
                        command="modalTemplate", settings=settings
                    )["points"]
                except RuntimeError as error:
                    if "Only " not in str(error):
                        raise
                    print(
                        json.dumps(dict(skipped=settings, reason=str(error))),
                        flush=True,
                    )
                    continue
                parameters = dict(renderer.initial if original is None else original)
                for i in range(32):
                    parameters[f"resolved_level_{i}"] = -72
                    parameters[f"resolved_turbulence_{i}"] = 1
                for i, point in enumerate(points):
                    parameters[f"resolved_frequency_{i}"] = point["frequency"]
                    parameters[f"resolved_level_{i}"] = -20
                for i in range(1, 7):
                    parameters[f"body_decay_active_{i}"] = 0
                # Neutral shared start, not detailed bars from the previous fit.
                if neutral:
                    parameters.update(
                        bloom_rate=3.3,
                        bloom_energy_acceleration=0.25,
                        body_brightness=-18,
                        body_excitation_centre=1000,
                        field_turbulence=0.65 * (1000 / 2500) ** 0.5,
                        field_turbulence_slope=0.5,
                        field_packet_spread=1.5,
                        field_phase_bandwidth=0.15,
                        field_wander_hz=0,
                        body_decay_seconds_0=7,
                        body_decay_seconds_7=1.2,
                    )
                yield settings, parameters


def shared_bounds():
    return dict(
        bloom_rate=(0.02, 16),
        bloom_energy_acceleration=(0, 1),
        body_brightness=(-48, 12),
        body_excitation_centre=(150, 4500),
        field_turbulence=(0.05, 3),
        field_turbulence_slope=(0, 1),
        field_packet_spread=(0.1, 6),
        field_phase_bandwidth=(0.001, 2),
        body_decay_seconds_0=(0.3, 20),
        body_decay_seconds_7=(0.15, 8),
    )


def run(args):
    if args.finalists < 1 or args.iterations < 1:
        raise ValueError("Positive finalist and iteration counts required")
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        mel = AuralossMel(reference, renderer.sample_rate)

        def save(directory, parameters, history):
            return checkpoint(
                renderer,
                shape,
                args.output / directory,
                f"{args.target.title()} - stretched series, coarse fit",
                parameters,
                reference,
                history,
            )

        def score(search):
            return mel.score(search.audio(search.parameters)) + 0.05 * np.linalg.norm(
                search.residual(search.parameters)
            )

        baseline = save(
            "baseline", renderer.initial, [dict(stage="current renderer baseline")]
        )
        print(json.dumps(dict(baseline_score=score(baseline))), flush=True)
        candidates = []
        original = (
            None
            if args.source is None
            else verify_candidate(renderer, args.source)["parameters"]
        )
        for index, (settings, parameters) in enumerate(
            layouts(
                renderer,
                args.roots,
                args.neutral,
                original,
                args.counts,
                args.stretches,
            )
        ):
            search = save(
                f"grid-{index}",
                parameters,
                [
                    dict(
                        stage="whole layout",
                        settings=settings,
                        shared_start="neutral" if args.neutral else "previous fit",
                        observation_knots=args.knots,
                    )
                ],
            )
            polish_coarse(search, args.knots)
            value = score(search)
            candidates.append((value, search))
            print(
                json.dumps(
                    dict(layout=settings, score=value, directory=str(search.output))
                ),
                flush=True,
            )
        if not candidates:
            raise ValueError("No complete modal layout fits the requested range")
        scales = {
            key: "log"
            for key in (
                "bloom_rate",
                "body_excitation_centre",
                "field_phase_bandwidth",
                "body_decay_seconds_0",
                "body_decay_seconds_7",
            )
        }
        finalists = []
        for rank, (_, source) in enumerate(
            sorted(candidates, key=lambda row: row[0])[: args.finalists]
        ):
            search = save(
                f"refined-{rank}",
                source.parameters,
                source.history + [dict(parent=str(source.output))],
            )
            for step in range(2):
                search.stage(
                    f"shared dynamics {step+1}",
                    shared_bounds(),
                    args.iterations,
                    difference_step=0.005,
                    parameter_scales=scales,
                )
                polish_coarse(search, args.knots)
            finalists.append((score(search), search))
        _, best = min(finalists, key=lambda row: row[0])
        result = save(
            "candidate", best.parameters, best.history + [dict(parent=str(best.output))]
        )
        seed = renderer.metadata["event"]["seed"]
        result.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
        perceptual_polish(result, mel, args.knots)
        verify_candidate(renderer, result.output)
        print(
            json.dumps(dict(final_score=score(result), directory=str(result.output))),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--roots", nargs="+", type=float, default=[137, 175, 205])
    parser.add_argument("--counts", nargs="+", type=int, default=[24, 32])
    parser.add_argument("--stretches", nargs="+", type=float, default=[0.25, 0.5, 0.75])
    parser.add_argument(
        "--source", type=Path, help="Validated checkpoint for shared-control warm start"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=36)
    parser.add_argument("--finalists", type=int, default=2)
    parser.add_argument(
        "--neutral",
        action="store_true",
        help="Reset shared dynamics instead of warm-starting them",
    )
    parser.add_argument(
        "--knots", nargs="+", type=float, default=[120, 400, 1500, 4000, 8000, 15000]
    )
    run(parser.parse_args())
