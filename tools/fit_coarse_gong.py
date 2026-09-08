"""Protected harmonic core, fixed mode centres, four broad observation knots.

No individual frequencies, mode levels, local noisiness, or local decay fitting.
The initial layout is chosen as a whole; subsequent changes are shared controls.
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


def layouts(renderer, original):
    for count in (16, 24, 32):
        for stretch in (0.4, 0.7, 1.0):
            settings = dict(
                family="harmonic",
                fundamental=121.8,
                count=count,
                stretch=stretch,
                harmonicCore=4,
                level=-20,
                rolloff=0,
            )
            try:
                points = renderer.request(command="modalTemplate", settings=settings)[
                    "points"
                ]
            except RuntimeError as error:
                if "Only " not in str(error):
                    raise
                print(json.dumps(dict(skipped=settings, reason=str(error))), flush=True)
                continue
            parameters = dict(original)
            for i in range(32):
                parameters[f"resolved_level_{i}"] = -72
                parameters[f"resolved_turbulence_{i}"] = 1
            for i, p in enumerate(points):
                parameters[f"resolved_frequency_{i}"] = p["frequency"]
                parameters[f"resolved_level_{i}"] = -20
            yield settings, parameters


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        mel = AuralossMel(reference, renderer.sample_rate)
        checkpoint(
            renderer,
            shape,
            args.output / "baseline",
            "Previous gong",
            renderer.initial,
            reference,
            [dict(stage="current renderer baseline")],
        )
        finalists = []
        for index, (settings, parameters) in enumerate(
            layouts(renderer, renderer.initial)
        ):
            search = checkpoint(
                renderer,
                shape,
                args.output / f"grid-{index}",
                "Gong harmonic core",
                parameters,
                reference,
                [dict(stage="whole layout", settings=settings)],
            )
            polish_coarse(search)
            score = mel.score(search.audio(search.parameters)) + 0.05 * np.linalg.norm(
                search.residual(search.parameters)
            )
            finalists.append((score, search))
            print(json.dumps(dict(layout=settings, score=score)), flush=True)
        finalists.sort(key=lambda item: item[0])
        bounds = dict(
            bloom_rate=(0.05, 16),
            bloom_energy_acceleration=(0, 1),
            body_brightness=(-60, 0),
            body_excitation_centre=(100, 3000),
            field_turbulence=(0.05, 2),
            field_turbulence_slope=(0.1, 1),
            field_packet_spread=(0.2, 5),
            field_phase_bandwidth=(0.0001, 1),
            body_decay_seconds_0=(1, 20),
            body_decay_seconds_7=(0.3, 10),
        )
        scales = {
            k: "log"
            for k in (
                "bloom_rate",
                "body_excitation_centre",
                "field_phase_bandwidth",
                "body_decay_seconds_0",
                "body_decay_seconds_7",
            )
        }
        fitted = []
        for rank, (_, source) in enumerate(finalists[:2]):
            search = checkpoint(
                renderer,
                shape,
                args.output / f"refined-{rank}",
                "Gong protected harmonic core",
                source.parameters,
                reference,
                source.history + [dict(parent=str(source.output))],
            )
            for stage in range(2):
                search.stage(
                    f"shared dynamics {stage+1}",
                    bounds,
                    24,
                    difference_step=0.005,
                    parameter_scales=scales,
                )
                polish_coarse(search)
            score = mel.score(search.audio(search.parameters)) + 0.05 * np.linalg.norm(
                search.residual(search.parameters)
            )
            fitted.append((score, search))
        _, chosen = min(fitted, key=lambda item: item[0])
        result = checkpoint(
            renderer,
            shape,
            args.output / "candidate",
            "Gong protected harmonic core",
            chosen.parameters,
            reference,
            chosen.history + [dict(parent=str(chosen.output))],
        )
        verify_candidate(renderer, result.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
