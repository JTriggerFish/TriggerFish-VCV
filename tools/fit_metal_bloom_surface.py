"""Fit exposed excitation, diffusion and packet controls to onset/bloom shape.

Frequencies, per-mode damping, gestures and output gain are fixed. Search the
nonlinearity as well as its coefficient: holding it at one is not justified.
"""

import argparse
import os
from pathlib import Path

import numpy as np
from scipy.stats import qmc
from scipy.optimize import least_squares

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomBasis
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.audio_io import AudioBuffer, write_wav


def polish(search):
    keys = [
        k
        for k, v in search.parameters.items()
        if k.startswith("resolved_level_") and v > -71.99
    ]
    basis = ObservationBasis(
        search.renderer, search.parameters, keys, search.seconds, search.seeds
    )
    cache = SpectralBloomBasis(basis, search.loss)
    low, high = 10 ** (-45 / 20), 10 ** (6 / 20)
    result = least_squares(
        lambda x: cache.evaluate(x)[0],
        np.clip(basis.amplitudes, low, high),
        jac=lambda x: cache.evaluate(x)[1],
        bounds=(low, high),
        max_nfev=200,
    )
    validation = cache.validate(basis, result.x)
    candidate = dict(
        search.parameters, **dict(zip(keys, (20 * np.log10(result.x)).tolist()))
    )
    before = float(np.linalg.norm(search.residual(search.parameters)))
    after = float(np.linalg.norm(search.residual(candidate)))
    if after < before:
        search.parameters = candidate
    search.history.append(
        dict(
            stage="exact STFT observation",
            before=before,
            after=after,
            selected=after < before,
            basis_error=validation,
            basis_validation=basis.validation,
        )
    )
    search.save()


def run(args):
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        seconds = 6
        loss = SpectralBloomLoss(
            aligned_reference(renderer, seconds), renderer.sample_rate
        )
        args.output.mkdir(parents=True, exist_ok=True)
        write_wav(
            args.output / "reference.wav",
            AudioBuffer(aligned_reference(renderer, seconds), renderer.sample_rate),
        )
        search = Search(renderer, loss, args.output, seconds, args.target + " bloom")
        bounds = dict(
            bloom_rate=(0.01, 16),
            bloom_energy_acceleration=(0, 1),
            bloom_energy_sensitivity=(0, 2),
            body_brightness=(-48, 0),
            body_excitation_centre=(100, 2500),
            field_turbulence=(0.05, 2),
            field_turbulence_slope=(0.2, 1),
            field_packet_spread=(0.2, 4),
            field_phase_bandwidth=(0.00001, 1),
            body_decay_seconds_0=(1, 30),
            body_decay_seconds_7=(0.3, 12),
        )
        initial = dict(search.parameters)
        if args.source:
            import json

            source = json.loads(args.source.read_text())
            search.parameters = source["parameters"]
            if set(search.parameters) != set(renderer.initial):
                raise ValueError(
                    "Warm start must contain the complete current parameter surface"
                )
            search.history.append(
                dict(
                    stage="explicit warm start",
                    source=str(args.source.resolve()),
                    parameters=search.parameters.copy(),
                    renderer_sha256=renderer.metadata["rendererSha256"],
                    previous_renderer=source.get("metadata", {}).get("rendererSha256"),
                    old_scores_reused=False,
                )
            )
        else:
            candidates = []
            for i, u in enumerate(qmc.LatinHypercube(9, seed=2431).random(96)):
                values = dict(
                    initial,
                    bloom_rate=float(0.03 * (16 / 0.03) ** u[0]),
                    bloom_energy_acceleration=float(u[1]),
                    bloom_energy_sensitivity=float(2 * u[8]),
                    body_brightness=float(-6 - 36 * u[2]),
                    body_excitation_centre=float(100 * 20 ** u[3]),
                    field_phase_bandwidth=float(0.0001 * 500 ** u[4]),
                    field_turbulence=float(0.05 + 1.2 * u[5]),
                    field_turbulence_slope=float(0.2 + 0.8 * u[6]),
                    field_packet_spread=float(0.2 + 3 * u[7]),
                )
                candidates.append((str(i), values))
            search.screen_candidates("joint excitation/diffusion starts", candidates)
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
        if args.decay_knot:
            if any(search.parameters[f"body_decay_active_{i}"] for i in range(1, 7)):
                raise ValueError(
                    "This final-shaping trial expects a two-endpoint source"
                )
            erb = lambda f: np.log1p(0.00437 * f)
            fraction = (erb(args.decay_knot) - erb(40)) / (erb(15000) - erb(40))
            seconds_at_knot = np.exp(
                (1 - fraction) * np.log(search.parameters["body_decay_seconds_0"])
                + fraction * np.log(search.parameters["body_decay_seconds_7"])
            )
            search.parameters.update(
                body_decay_active_1=1,
                body_decay_frequency_1=args.decay_knot,
                body_decay_seconds_1=float(seconds_at_knot),
            )
            bounds = {
                k: bounds[k] for k in ("body_decay_seconds_0", "body_decay_seconds_7")
            }
            bounds["body_decay_seconds_1"] = (0.3, 20)
            scales["body_decay_seconds_1"] = "log"
            search.history.append(
                dict(
                    stage="explicit final shared damping knot",
                    frequency_hz=args.decay_knot,
                    initial_seconds=float(seconds_at_knot),
                )
            )
        if args.dynamics_only:
            bounds = {
                k: v
                for k, v in bounds.items()
                if k
                in (
                    "bloom_rate",
                    "bloom_energy_acceleration",
                    "body_brightness",
                    "body_excitation_centre",
                    "body_decay_seconds_0",
                    "body_decay_seconds_7",
                )
            }
        search.stage(
            "onset and bloom",
            bounds,
            32,
            difference_step=0.005,
            parameter_scales=scales,
        )
        # Positive observation levels are independent of the state evolution.
        polish(search)
        search.stage(
            "joint shape refinement",
            bounds,
            30,
            difference_step=0.003,
            parameter_scales=scales,
        )
        polish(search)
        search.save()
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--dynamics-only", action="store_true")
    parser.add_argument("--decay-knot", type=float)
    args = parser.parse_args()
    if args.decay_knot is not None and not 40 < args.decay_knot < 15000:
        parser.error("Interior damping knot must lie between 40 Hz and 15 kHz")
    run(args)
