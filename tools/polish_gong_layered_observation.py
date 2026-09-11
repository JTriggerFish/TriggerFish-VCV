"""Fit a four-point observation curve using validated exact-render STFT bases.

Centres and dynamics are fixed. Solver coordinates become explicit painted
levels; the synthesizer has no hidden curve or extra processing stage.
"""

import argparse
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch
from scipy.optimize import minimize

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomBasis
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_gong_layered_bloom import LayeredBloomLoss


def paint(base, gains):
    """Smooth dB interpolation over four frequency anchors, not individual fits."""
    centres = np.log([240, 600, 3000, 12000])
    result = dict(base)
    for i in range(32):
        if base[f"resolved_level_{i}"] <= -71.99:
            continue
        f = np.log(base[f"resolved_frequency_{i}"])
        j = np.clip(np.searchsorted(centres, f) - 1, 0, 2)
        x = np.clip((f - centres[j]) / (centres[j + 1] - centres[j]), 0, 1)
        x = x * x * (3 - 2 * x)
        gain = (1 - x) * gains[j] + x * gains[j + 1]
        result[f"resolved_level_{i}"] = float(
            np.clip(base[f"resolved_level_{i}"] + gain, -60, 6)
        )
    return result


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        user = json.loads((args.source / "baseline/search.json").read_text())[
            "parameters"
        ]
        ref = aligned_reference(r, 6)
        loss = LayeredBloomLoss(ref, r.render(user, 6), r.sample_rate)
        rows = json.loads((args.source / "screen.json").read_text())
        dynamic = (
            json.loads(args.parameters.read_text())["parameters"]
            if args.parameters
            else rows[args.trial]["parameters"]
        )
        # Restore original bars before constructing the independent basis.
        p = dict(dynamic)
        if args.body_observation is not None:
            p["field_gain"] = args.body_observation
        keys = [f"resolved_level_{i}" for i in range(32)]
        p.update({k: user[k] for k in keys})
        basis = ObservationBasis(r, p, keys, 6, (1675, 1982))
        proxy = SimpleNamespace(
            rate=r.sample_rate,
            edges=[80, 300, 900, 3000, 7000, 14000],
            regions=list(zip(loss.times[:-1], loss.times[1:])),
            active=np.ones(5, dtype=bool),
        )
        cache = SpectralBloomBasis(basis, proxy)

        def diagnostics(gains):
            params = paint(p, gains)
            vector = np.r_[1.0, 10 ** (np.array([params[k] for k in keys]) / 20)]
            power = (cache.matrices @ vector) @ vector
            db = 10 * np.log10(np.maximum(power, 1e-10))
            error = db - loss.target
            rms = lambda x: np.sqrt(np.mean(x * x, axis=tuple(range(1, x.ndim))))
            values = 2 * rms(error[:, 0]) + rms(error[:, 1:3]) + rms(error[:, 3:])
            values += rms(error[:, 3:, 3:8] - error[:, 3:, 1:2])
            return float(values.mean()), params, db

        offset = 20 * np.log10(p["field_gain"] / user["field_gain"])
        initial = [-16 - offset, -offset, 16 - offset, 36 - offset]
        result = minimize(
            lambda x: diagnostics(x)[0],
            initial,
            method="Powell",
            bounds=[(-30, -10), (-16, 12), (-8, 32), (8, 48)],
            options=dict(maxfev=500, xtol=0.005, ftol=0.0001),
        )
        score, parameters, predicted = diagnostics(result.x)
        # Reject an incorrect affine/cache shortcut before saving anything.
        actual = np.array(
            [loss.envelopes(r.render(parameters, 6, s)) for s in (1675, 1982)]
        )
        error = float(abs(actual - predicted).max())
        if error > 0.02:
            raise ValueError(f"Exact band-envelope verification failed: {error} dB")
        search = checkpoint(
            r,
            loss,
            args.output,
            "Gong — tuned body and bloom",
            parameters,
            ref,
            [
                dict(
                    stage="four-point grouped observation",
                    dynamic_trial=args.trial,
                    dynamic_checkpoint=(
                        str(args.parameters) if args.parameters else None
                    ),
                    body_observation=p["field_gain"],
                    dynamics={
                        k: p[k]
                        for k in (
                            "bloom_rate",
                            "bloom_energy_acceleration",
                            "bloom_energy_sensitivity",
                        )
                    },
                    curve_hz=[240, 600, 3000, 12000],
                    curve_db=result.x.tolist(),
                    solver="Powell on exact STFT cross-power",
                    evaluations=int(result.nfev),
                    basis_validation=basis.validation,
                    cache_error_db=error,
                    objective=loss.specification,
                    score=score,
                )
            ],
        )
        verify_candidate(r, args.output)
        audit = [
            dict(
                seed=s,
                baseline=loss.diagnostics(r.render(user, 6, s)),
                candidate=loss.diagnostics(r.render(parameters, 6, s)),
            )
            for s in (1675, 1982, 2586, 3276)
        ]
        (args.output / "audit.json").write_text(json.dumps(audit, indent=2))
        print(
            json.dumps(
                dict(score=score, curve=result.x.tolist(), error=error, audit=audit)
            ),
            flush=True,
        )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--trial", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parameters", type=Path)
    parser.add_argument("--body-observation", type=float)
    run(parser.parse_args())
