"""Optional actual-Wasm contact observation fitting; never publishes presets."""

import argparse
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from scipy.optimize import minimize

from review_instrument_fit import band_audit
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.observation_energy_gate import EARLY_BINS
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search

BOUNDS = {
    "direct_gain": (0.5, 2.0),
}


def early_power_db(audio, rate):
    """Fixed, unnormalized early RMS bins; reject malformed rendered audio."""
    if (
        audio.ndim != 1
        or len(audio) < round(rate * 0.1)
        or not np.isfinite(audio).all()
    ):
        raise ValueError("Expected finite mono audio covering all attack bins")
    values = np.array(
        [
            10
            * np.log10(
                max(
                    float(np.mean(audio[round(a * rate) : round(b * rate)] ** 2)), 1e-20
                )
            )
            for a, b in EARLY_BINS
        ]
    )
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite attack power")
    return values


def derivative(values, function, step=0.005):
    """Central finite differences in bounded normalized control coordinates."""
    columns = []
    for index in range(len(values)):
        minus, plus = values.copy(), values.copy()
        minus[index] = max(0, values[index] - step)
        plus[index] = min(1, values[index] + step)
        columns.append(
            (function(plus) - function(minus)) / (plus[index] - minus[index])
        )
    return np.asarray(columns).T


def fit_contact(renderer, parameters, reference, iterations, tolerance):
    """Constrain every short RMS and first-30-ms band error for each seed."""
    keys = list(BOUNDS)
    low, high = np.array(list(BOUNDS.values())).T
    descriptors = {row["key"]: row for row in renderer.metadata["descriptors"]}
    for key, minimum, maximum in zip(keys, low, high):
        if (
            descriptors[key]["minimum"] > minimum
            or descriptors[key]["maximum"] < maximum
        ):
            raise ValueError(f"Contact fitting bounds exceed exposed range: {key}")
    seeds = (None, (renderer.metadata["event"]["seed"] + 101) & 0xFFFFFFFF)
    target = early_power_db(reference, renderer.sample_rate)
    cache, best = {}, {}

    def unpack(values):
        return dict(
            parameters, **dict(zip(keys, (low + values * (high - low)).tolist()))
        )

    def errors(values):
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite trial controls")
        token = values.tobytes()
        if token not in cache:
            rows = []
            for seed in seeds:
                audio = renderer.render(unpack(values), 0.2, seed)
                rows.extend(early_power_db(audio, renderer.sample_rate) - target)
                rows.extend(
                    band_audit(reference, audio, renderer.sample_rate)["regions"][0][
                        "difference_db"
                    ]
                )
            cache[token] = np.asarray(rows)
            if not np.isfinite(cache[token]).all():
                raise ValueError("Nonfinite band errors")
        return cache[token]

    def objective(values):
        error = errors(values)
        score = float(np.mean(error**2))
        if np.max(np.abs(error)) <= tolerance and score < best.get("score", np.inf):
            best.update(
                score=score, parameters=unpack(values), errors_db=error.tolist()
            )
        return score

    def constraints(values):
        error = errors(values)
        return np.r_[tolerance - error, tolerance + error]

    # Explicit diagnostic start, not an undisclosed clamp of the input patch.
    start = (np.array([1.4, 1000, 4000, 3]) - low) / (high - low)
    result = minimize(
        objective,
        start,
        jac=lambda x: derivative(x, objective),
        bounds=[(0, 1)] * 4,
        method="SLSQP",
        constraints=[
            dict(type="ineq", fun=constraints, jac=lambda x: derivative(x, constraints))
        ],
        options=dict(maxiter=iterations, ftol=1e-6),
    )
    if not np.isfinite(result.x).all() or not np.isfinite(result.fun):
        raise ValueError("Nonfinite optimizer result")
    if not best:
        raise RuntimeError("No feasible spectral/energy attack; no fit saved")
    return best, dict(
        stage="contact observation spectral and early energy constraints",
        algorithm="SLSQP",
        bounds=BOUNDS,
        initial_controls=dict(zip(keys, [1.4, 1000, 4000, 3])),
        difference_step_fraction=0.005,
        tolerance_db=tolerance,
        training_seeds=seeds,
        maximum_iterations=iterations,
        iterations=int(result.nit),
        ftol=1e-6,
        solver_success=bool(result.success),
        solver_message=str(result.message),
        evaluation_seconds=0.2,
        early_bins_seconds=EARLY_BINS,
        band_region_seconds=[0, 0.03],
        bands_hz=[[40, 250], [250, 1000], [1000, 4000], [4000, 16000]],
        fixed_parameters={k: v for k, v in parameters.items() if k not in keys},
    )


def refine(args):
    if (
        args.start.resolve() == args.output.resolve()
        or (args.output / "search.json").exists()
    ):
        raise ValueError("Output must be separate and contain no existing search.json")
    if (
        args.iterations < 1
        or not np.isfinite(args.tolerance_db)
        or args.tolerance_db <= 0
    ):
        raise ValueError("Iterations and finite tolerance must be positive")
    if not np.isfinite(args.seconds) or args.seconds < 3.1:
        raise ValueError(
            "Full validation duration must be finite and at least 3.1 seconds"
        )
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, args.start)
        short = aligned_reference(renderer, 0.2)
        best, record = fit_contact(
            renderer, saved["parameters"], short, args.iterations, args.tolerance_db
        )
        reference = aligned_reference(renderer, args.seconds)
        loss = MetallicBalanceLoss(
            reference, renderer.sample_rate, contrast_weighting="erb", fast_attack=True
        )
        for seed in record["training_seeds"]:
            audio = renderer.render(best["parameters"], args.seconds, seed)
            errors = np.r_[
                early_power_db(audio, renderer.sample_rate)
                - early_power_db(reference, renderer.sample_rate),
                band_audit(reference, audio, renderer.sample_rate)["regions"][0][
                    "difference_db"
                ],
            ]
            if (
                not np.isfinite(errors).all()
                or np.max(abs(errors)) > args.tolerance_db + 0.02
            ):
                raise ValueError("Fresh full render failed contact constraints")
        args.output.mkdir(parents=True, exist_ok=True)
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        search = Search(
            renderer,
            loss,
            args.output,
            args.seconds,
            args.target.title(),
            tuple(record["training_seeds"]),
        )
        search.parameters = best["parameters"]
        search.history = [
            dict(record, best_score=best["score"], actual_errors_db=best["errors_db"])
        ]
        search.save()
        verify_candidate(renderer, args.output)
        print(json.dumps(search.history[-1]), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("ride", "crash", "gong", "hihat"))
    parser.add_argument("start", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seconds", type=float, default=12)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--tolerance-db", type=float, default=3.5)
    refine(parser.parse_args())
