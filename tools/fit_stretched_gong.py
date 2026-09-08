"""Start a gong fit from the UI's stretched harmonic generator, not old ridges.

Screen 24 explicit grids, polish the best three with exact-render observation
bases, then refine shared dynamics. Archive the existing preset as the baseline.
No gesture/gain normalization, added damping knots or per-mode damping fits.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from fit_metal_bloom_surface import polish


def seed_parameters(renderer, original, settings):
    # Calling the same JS generator prevents another interpretation of stretch.
    points = renderer.request(
        command="modalTemplate",
        settings=dict(
            settings, family="harmonic", minimumFrequency=1, maximumFrequency=15000
        ),
    )["points"]
    old = sorted(
        (original[f"resolved_frequency_{i}"], original[f"resolved_level_{i}"])
        for i in range(32)
        if original[f"resolved_level_{i}"] > -71.99
    )
    frequency, level = np.array(old).T
    result = dict(original)
    for i in range(32):
        result[f"resolved_level_{i}"] = -72
        result[f"resolved_turbulence_{i}"] = 1
    for i, point in enumerate(points):
        result[f"resolved_frequency_{i}"] = point["frequency"]
        result[f"resolved_level_{i}"] = float(
            np.interp(np.log(point["frequency"]), np.log(frequency), level)
        )
    return result


def checkpoint(renderer, loss, directory, name, parameters, reference, history):
    directory.mkdir(parents=True, exist_ok=True)
    search = Search(renderer, loss, directory, 6, name)
    search.parameters = dict(parameters)
    search.history = list(history)
    search.save()
    write_wav(directory / "reference.wav", AudioBuffer(reference, renderer.sample_rate))
    return search


def screen(renderer, original, shape, mel, output):
    candidates = []
    for root in (55.0, 82.4069, 110.0, 123.4708):
        for count in (16, 24, 32):
            for top in (10000.0, 14000.0):
                try:
                    stretch = renderer.request(
                        command="modalTemplateStretch",
                        settings=dict(
                            fundamental=root,
                            count=count,
                            topFrequency=top,
                            harmonicCore=4,
                        ),
                    )["stretch"]
                except RuntimeError as error:
                    if "Top frequency is outside" not in str(error):
                        raise
                    print(
                        json.dumps(
                            dict(
                                skipped=dict(root=root, count=count, top=top),
                                reason=str(error),
                            )
                        ),
                        flush=True,
                    )
                    continue
                settings = dict(
                    fundamental=root,
                    count=count,
                    stretch=stretch,
                    harmonicCore=4,
                )
                parameters = seed_parameters(renderer, original, settings)
                audio = renderer.render(parameters, 6)
                envelope = float(np.linalg.norm(shape.residual(audio)))
                perceptual = mel.score(audio)
                row = dict(
                    settings=settings,
                    parameters=parameters,
                    bloom=envelope,
                    mel=perceptual,
                    rank_score=perceptual + 0.05 * envelope,
                )
                candidates.append(row)
                print(
                    json.dumps({k: v for k, v in row.items() if k != "parameters"}),
                    flush=True,
                )
    candidates.sort(key=lambda row: row["rank_score"])
    (output / "grid-screen.json").write_text(json.dumps(candidates, indent=2))
    return candidates[:3]


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
            "Gong prior preset",
            renderer.initial,
            reference,
            [dict(stage="current-renderer baseline")],
        )
        candidates = screen(renderer, renderer.initial, shape, mel, args.output)
        finalists = []
        for i, row in enumerate(candidates):
            search = checkpoint(
                renderer,
                shape,
                args.output / f"grid-{i}",
                "Gong stretched harmonic start",
                row["parameters"],
                reference,
                [
                    dict(
                        stage="UI generator seed",
                        settings=row["settings"],
                        observation_start="log-frequency interpolation of prior observation bars",
                    )
                ],
            )
            polish(search)
            samples = search.audio(search.parameters)
            score = mel.score(samples) + 0.05 * np.linalg.norm(shape.residual(samples))
            finalists.append((score, search))
        _, chosen = min(finalists, key=lambda item: item[0])
        search = checkpoint(
            renderer,
            shape,
            args.output / "candidate",
            "Gong stretched harmonic fit",
            chosen.parameters,
            reference,
            chosen.history + [dict(selected_from=str(chosen.output))],
        )
        bounds = dict(
            bloom_rate=(0.01, 16),
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
            key: "log"
            for key in (
                "bloom_rate",
                "body_excitation_centre",
                "field_phase_bandwidth",
                "body_decay_seconds_0",
                "body_decay_seconds_7",
            )
        }
        for name in ("stretched-grid dynamics", "stretched-grid refinement"):
            search.stage(
                name, bounds, 24, difference_step=0.005, parameter_scales=scales
            )
            polish(search)
        search.save()
        verify_candidate(renderer, search.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
