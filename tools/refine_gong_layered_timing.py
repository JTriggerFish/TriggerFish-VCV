"""Screen transport, then fit smooth observation balances with exact WASM bases.

No modal frequencies, damping multipliers, gestures or output gains are fitted.
Seven broad observation coordinates are saved as ordinary explicit painted bars.
Shape-only screening is not acceptance; final scoring uses reference levels too.
"""

import argparse
import json
import os
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch
from scipy.optimize import minimize

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomBasis
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint

ANCHORS = np.array([120, 360, 900, 2000, 4500, 9000, 14000])


def paint(base, gains, absolute=False, anchors=ANCHORS):
    """Broad smooth log-frequency balance; never change active mode membership."""
    p = dict(base)
    for i in range(32):
        key = f"resolved_level_{i}"
        if p[key] <= -71.99:
            continue
        f = np.log(p[f"resolved_frequency_{i}"])
        j = int(np.clip(np.searchsorted(np.log(anchors), f) - 1, 0, len(anchors) - 2))
        x = np.clip(
            (f - np.log(anchors[j])) / np.log(anchors[j + 1] / anchors[j]), 0, 1
        )
        x = x * x * (3 - 2 * x)
        level = 0 if absolute else p[key]
        p[key] = float(np.clip(level + (1 - x) * gains[j] + x * gains[j + 1], -60, 6))
    return p


def screen(r, loss, output, local=False, decay=False):
    rows = []
    trials = [
        (
            r.initial["bloom_rate"],
            r.initial["bloom_energy_acceleration"],
            r.initial["bloom_energy_sensitivity"],
        )
    ]
    trials += list(
        product((3, 4, 5, 6), (0.04, 0.08, 0.12, 0.16), (0.2, 0.3, 0.4))
        if local
        else product((4, 8, 12, 16), (0.02, 0.1, 0.3, 0.6, 1), (0, 0.3, 0.7))
    )
    trials = [
        (*trial, r.initial["body_decay_seconds_0"], r.initial["body_decay_seconds_7"])
        for trial in trials
    ]
    if decay:
        trials = [trials[0]] + list(
            product((4, 6), (0.08, 0.15), (0.3,), (3, 5, 7, 10), (1.5, 2.14, 3))
        )
    for index, (rate, concentration, sensitivity, low_decay, high_decay) in enumerate(
        trials
    ):
        p = dict(
            r.initial,
            bloom_rate=rate,
            bloom_energy_acceleration=concentration,
            bloom_energy_sensitivity=sensitivity,
            body_decay_seconds_0=low_decay,
            body_decay_seconds_7=high_decay,
        )
        db = loss.envelopes(r.render(p, 6, 1675))
        rows.append(
            dict(
                parameters=p,
                shape=loss.score_db(db, True),
                absolute=loss.score_db(db),
                envelopes_db=db.tolist(),
            )
        )
        if index % 10 == 0:
            print(
                json.dumps(dict(screen=index, best=min(x["shape"] for x in rows))),
                flush=True,
            )
    rows.sort(key=lambda row: row["shape"])
    (output / "screen.json").write_text(json.dumps(rows, indent=2))
    return rows


def polish(
    r,
    loss,
    ref,
    row,
    output,
    absolute=False,
    reuse_cache=None,
    anchors=ANCHORS,
    fit_seeds=1,
):
    anchors = np.asarray(anchors, dtype=float)
    if (
        not 3 <= len(anchors) <= 12
        or not np.isfinite(anchors).all()
        or anchors[0] <= 0
        or not np.all(np.diff(anchors) > 0)
    ):
        raise ValueError("Expected 3–12 finite, positive, increasing curve anchors")
    p = row["parameters"]
    keys = [
        f"resolved_level_{i}" for i in range(32) if p[f"resolved_level_{i}"] > -71.99
    ]
    if reuse_cache:
        saved = verify_candidate(r, reuse_cache)
        if any(v != saved["parameters"][k] for k, v in p.items() if k not in keys):
            raise ValueError("Cached observation dynamics differ from the source")
        with np.load(reuse_cache / "observation-cache.npz") as archive:
            matrices = archive["matrices"].copy()
        cache = SimpleNamespace(matrices=matrices)
        basis = SimpleNamespace(
            validation=[dict(reused_from=str(reuse_cache))],
            amplitudes=10 ** (np.array([p[k] for k in keys]) / 20),
        )
    else:
        basis = ObservationBasis(r, p, keys, 6, (1675, 1982))
        cache = SpectralBloomBasis(basis, loss)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "observation-cache.npz",
        matrices=cache.matrices,
        target=loss.target,
        lower=loss.lower,
        baseline_amplitudes=basis.amplitudes,
    )

    def evaluate(gains):
        params = paint(p, gains, absolute, anchors)
        v = np.r_[1.0, 10 ** (np.array([params[k] for k in keys]) / 20)]
        db = 10 * np.log10(np.maximum(cache.matrices @ v @ v, loss.floor))
        score = loss.score_db(db[:fit_seeds])
        # Weak broad-shape regularization, never independent centre placement.
        regularizer = 0.0005 * np.mean(np.diff(gains, 2) ** 2)
        return score + regularizer, params, db

    best = None
    starts = (np.zeros(len(ANCHORS)), np.array([-8, 4, -4, -8, -6, -4, -4]))
    if absolute:
        starts = (
            np.full(len(ANCHORS), -12.0),
            np.array([-35, -12, -12, -6, -6, -6, 3.0]),
        )
    for original_start in starts:
        start = np.interp(np.log(anchors), np.log(ANCHORS), original_start)
        result = minimize(
            lambda x: evaluate(x)[0],
            start,
            method="Powell",
            bounds=[(-60, 6) if absolute else (-30, 20)] * len(anchors),
            options=dict(maxfev=1800, xtol=0.02, ftol=1e-5),
        )
        if best is None or result.fun < best.fun:
            best = result
    score, params, predicted = evaluate(best.x)
    actual = np.array([loss.envelopes(r.render(params, 6, s)) for s in (1675, 1982)])
    cache_error = float(np.max(abs(actual - predicted)))
    if not np.isfinite(cache_error) or cache_error > 0.02:
        raise ValueError(f"Actual-render cache verification failed: {cache_error} dB")
    history = [
        dict(
            stage="transport plus smooth observation",
            parameters=params,
            anchors_hz=list(anchors),
            gains_db=best.x.tolist(),
            absolute_observation_curve=absolute,
            fitting_seeds=[1675, 1982][:fit_seeds],
            cache_error_db=cache_error,
            basis_validation=basis.validation,
            specification=loss.specification,
            score=score,
        )
    ]
    checkpoint(
        r, loss, output, "Gong — separated body and upper bloom", params, ref, history
    )
    verify_candidate(r, output)
    summary = dict(
        score=score,
        parameters=params,
        gains_db=best.x.tolist(),
        cache_error_db=cache_error,
        directory=str(output),
    )
    print(
        json.dumps({k: v for k, v in summary.items() if k != "parameters"}), flush=True
    )
    return summary


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        if args.source:
            r.initial = json.loads(args.source.read_text())["parameters"]
        ref = aligned_reference(r, 6)
        loss = LayeredBandLoss(ref, r.sample_rate, args.audibility)
        if args.polish_source:
            rows = [dict(parameters=r.initial)]
        elif not args.polish_only:
            checkpoint(
                r, loss, args.output / "baseline", "Incoming gong", r.initial, ref, []
            )
            rows = screen(r, loss, args.output, args.local, args.decay)
        else:
            rows = json.loads((args.output / "screen.json").read_text())
        results = []
        for index in args.ranks:
            results.append(
                polish(
                    r,
                    loss,
                    ref,
                    rows[index],
                    args.output / f"trial-{index}",
                    args.absolute,
                    args.reuse_cache,
                    args.anchors,
                    args.fit_seeds,
                )
            )
        (args.output / "polished.json").write_text(json.dumps(results, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--polish-only", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--local", action="store_true")
    parser.add_argument("--decay", action="store_true")
    parser.add_argument("--absolute", action="store_true")
    parser.add_argument("--polish-source", action="store_true")
    parser.add_argument("--audibility", action="store_true")
    parser.add_argument("--reuse-cache", type=Path)
    parser.add_argument("--anchors", type=float, nargs="+", default=ANCHORS.tolist())
    parser.add_argument("--fit-seeds", type=int, choices=[1, 2], default=1)
    parser.add_argument("--ranks", type=int, nargs="*", default=[0, 1, 2])
    run(parser.parse_args())
