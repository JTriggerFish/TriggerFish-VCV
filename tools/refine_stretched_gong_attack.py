"""Test missing measured low ridges after a stretched-grid gong initialization.

Keep the handle budget. Retune three low centres and exchange two nearly
inaudible handles for measured 247/375-Hz attack ridges. Then refit observation
and shared dynamics before comparing the alternative against the strict grid.
No hidden resonators, normalization or per-mode damping changes.
"""

import argparse
import os
from pathlib import Path

import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from fit_metal_bloom_surface import polish


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer,
            SpectralBloomLoss(reference, renderer.sample_rate),
            args.output,
            6,
            "Gong - stretched grid with measured low ridges",
        )
        search.parameters = dict(saved["parameters"])
        search.history.append(dict(parent=str(args.source.resolve())))
        quiet = [
            i
            for i in range(32)
            if search.parameters[f"resolved_frequency_{i}"] > 800
            and -71.99 < search.parameters[f"resolved_level_{i}"] < -40
        ]
        quiet.sort(key=lambda i: search.parameters[f"resolved_level_{i}"])
        if len(quiet) < 2:
            raise ValueError(
                "This hypothesis needs two quiet upper handles to exchange"
            )
        # Measured first-400-ms Welch ridges. Do not let a broad Mel objective
        # trade their pitch for a lower aggregate error. The 247-Hz ridge is
        # clearer later in the recording than the initial 344/375-Hz pair.
        centres = {0: 121.8, 1: 343.9, 2: 553.8, quiet[0]: 247.0, quiet[1]: 374.1}
        for index, frequency in centres.items():
            search.parameters[f"resolved_frequency_{index}"] = frequency
        for index in quiet[:2]:
            search.parameters[f"resolved_level_{index}"] = -25
        search.history.append(
            dict(
                stage="measured low ridge hypothesis",
                centres=centres,
                method="first 400 ms Welch: 16384 Hann samples, 65536 FFT; manual ridge selection",
                untouched_seed_parameters=saved["parameters"],
            )
        )
        polish(search)
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
        search.stage(
            "shared dynamics after ridge placement",
            bounds,
            24,
            difference_step=0.005,
            parameter_scales=scales,
        )
        polish(search)
        search.save()
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
