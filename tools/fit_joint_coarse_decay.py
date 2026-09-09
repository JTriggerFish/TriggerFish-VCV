"""Joint broad prominence/damping fit: avoid a level-versus-decay coordinate trap.

Fixed modal centres, local noisiness, source, gesture and output gains. Only six
broad positive observation amplitudes and a sparse shared T60 curve are free.
Every residual and finite difference uses the exact Wasm, no surrogate voice.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
import torch

from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_coarse_metal import KNOTS, damping_variant


def fit(search, frequencies, iterations):
    if iterations < 1:
        raise ValueError("Positive iteration budget required")
    original = dict(search.parameters)
    prepared = damping_variant(original, frequencies)
    indices = [i for i in range(32) if original[f"resolved_level_{i}"] > -71.99]
    weights = interpolation_weights(
        [original[f"resolved_frequency_{i}"] for i in indices], KNOTS
    )
    levels = np.array([original[f"resolved_level_{i}"] for i in indices])
    amplitudes = np.linalg.lstsq(weights, 10 ** (levels / 20), rcond=None)[0]
    decay_keys = [
        "body_decay_seconds_0",
        *[f"body_decay_seconds_{i}" for i in range(1, len(frequencies) + 1)],
        "body_decay_seconds_7",
    ]
    low = np.r_[np.full(6, 10 ** (-45 / 20)), np.full(len(decay_keys), np.log(0.1))]
    high = np.r_[np.full(6, 10 ** (6 / 20)), np.full(len(decay_keys), np.log(30))]
    start = np.clip(
        np.r_[amplitudes, np.log([prepared[k] for k in decay_keys])], low, high
    )

    def parameters(x):
        result = dict(prepared)
        for i, amplitude in zip(indices, weights @ x[:6]):
            result[f"resolved_level_{i}"] = float(
                np.clip(20 * np.log10(amplitude), -45, 6)
            )
        result.update(dict(zip(decay_keys, np.clip(np.exp(x[6:]), 0.1, 30).tolist())))
        return result

    def residual(x):
        return search.residual(parameters(x))

    counter = 0

    def jacobian(x):
        nonlocal counter
        counter += 1
        columns = []
        for i in range(len(x)):
            step = 0.005 if i < 6 else 0.02
            a, b = x.copy(), x.copy()
            a[i], b[i] = max(low[i], x[i] - step), min(high[i], x[i] + step)
            columns.append((residual(b) - residual(a)) / (b[i] - a[i]))
        record = dict(
            iteration=counter,
            score=float(np.linalg.norm(residual(x))),
            parameters=parameters(x),
        )
        (search.output / "joint-progress.json").write_text(
            json.dumps(record), encoding="utf8"
        )
        print(
            json.dumps({k: v for k, v in record.items() if k != "parameters"}),
            flush=True,
        )
        return np.array(columns).T

    before = float(np.linalg.norm(search.residual(original)))
    result = least_squares(
        residual,
        start,
        jac=jacobian,
        bounds=(low, high),
        x_scale="jac",
        max_nfev=iterations,
        ftol=0.0005,
        xtol=0.0005,
        gtol=0.001,
    )
    after = float(np.linalg.norm(residual(result.x)))
    accepted = np.isfinite(after) and after < before
    search.parameters = parameters(result.x) if accepted else original
    search.history.append(
        dict(
            stage="joint broad prominence and shared damping",
            knots_hz=KNOTS,
            decay_frequencies=frequencies,
            amplitude_bounds_db=[-45, 6],
            decay_bounds_seconds=[0.1, 30],
            amplitude_probe=0.005,
            log_decay_probe=0.02,
            before=before,
            trial_score=after,
            after=after if accepted else before,
            accepted=bool(accepted),
            solver=str(result.message),
            evaluations=int(result.nfev),
        )
    )
    search.save()


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        search = checkpoint(
            renderer,
            shape,
            args.output,
            f"{args.target.title()} - joint coarse prominence and decay",
            source["parameters"],
            reference,
            source["history"] + [dict(parent=str(args.source))],
        )
        fit(search, args.frequencies, args.iterations)
        seed = renderer.metadata["event"]["seed"]
        search.seeds = (seed, (seed + 101) & 0xFFFFFFFF)
        polish_coarse_mel(search, AuralossMel(reference, renderer.sample_rate), KNOTS)
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frequencies", type=float, nargs="+", default=[400, 800])
    parser.add_argument("--iterations", type=int, default=24)
    run(parser.parse_args())
