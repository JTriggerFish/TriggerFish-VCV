"""Optional actual-Wasm observation fitting with explicit early-energy gates.

Example after dev.ps1 build-workbench and the SDK environment is initialized:
  python tools/refine_workbench_attack_gate.py ride START_DIRECTORY OUTPUT_DIRECTORY
This never edits presets, serves pages, normalizes audio, or changes strike inputs.
"""

import argparse
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from scipy.optimize import minimize
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.observation_energy_gate import ObservationEnergyGate
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.torch_metallic_loss import TorchMetallicLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def gated_weights(basis, measurement, gate, iterations):
    """SLSQP, retaining the best feasible evaluation even if the last one fails."""
    lower, upper = 10 ** (-45 / 20), 10 ** (6 / 20)
    if iterations < 1:
        raise ValueError("Iterations must be positive")
    if (
        not np.isfinite(basis.amplitudes).all()
        or np.any(basis.amplitudes < lower)
        or np.any(basis.amplitudes > upper)
    ):
        raise ValueError(
            "Active input levels must lie within the explicit -45 to +6 dB bounds"
        )
    prepared = [(torch.tensor(a), torch.tensor(c)) for a, c in basis.bases.values()]
    initial = torch.tensor(basis.amplitudes)
    best = {"score": float("inf"), "weights": None}
    cache = {}

    def objective(values):
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite trial amplitudes")
        key = values.tobytes()
        if cache.get("key") == key:
            return cache["result"]
        weights = torch.tensor(values, requires_grad=True)
        total = 0.0
        for audio, columns in prepared:
            value = measurement(audio + (weights - initial) @ columns) / len(prepared)
            value.backward()
            total += float(value.detach())
        gradient = weights.grad.numpy().copy()
        if not np.isfinite(total) or not np.isfinite(gradient).all():
            raise ValueError("Nonfinite objective or gradient")
        if gate.evaluate(values)[0].min() >= -1e-6 and total < best["score"]:
            best.update(score=total, weights=values.copy())
        cache.update(key=key, result=(total, gradient))
        return total, gradient

    result = minimize(
        objective,
        basis.amplitudes,
        jac=True,
        method="SLSQP",
        bounds=[(lower, upper)] * len(basis.keys),
        constraints=[
            dict(
                type="ineq",
                fun=lambda x: gate.evaluate(x)[0],
                jac=lambda x: gate.evaluate(x)[1],
            )
        ],
        options=dict(maxiter=iterations, ftol=1e-6),
    )
    if not np.isfinite(result.x).all() or not np.isfinite(result.fun):
        raise ValueError("Nonfinite optimizer result")
    if best["weights"] is None:
        raise RuntimeError("No feasible early-energy candidate; no output fit saved")
    return best["weights"], dict(
        algorithm="SLSQP",
        solver_message=str(result.message),
        solver_success=bool(result.success),
        iterations=int(result.nit),
        best_feasible_score=float(np.sqrt(best["score"])),
    )


def refine(arguments):
    if (
        arguments.start.resolve() == arguments.output.resolve()
        or (arguments.output / "search.json").exists()
    ):
        raise ValueError(
            "Output must be separate from the input and contain no existing search.json"
        )
    if not np.isfinite(arguments.seconds) or arguments.seconds <= 0:
        raise ValueError("Render duration must be finite and positive")
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], arguments.target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, arguments.start)
        reference = aligned_reference(renderer, arguments.seconds)
        loss = MetallicBalanceLoss(
            reference, renderer.sample_rate, contrast_weighting="erb", fast_attack=True
        )
        parameters = saved["parameters"]
        keys = [
            k
            for k, v in parameters.items()
            if k.startswith("resolved_level_") and v > -71.99
        ]
        seeds = (None, (renderer.metadata["event"]["seed"] + 101) & 0xFFFFFFFF)
        basis = ObservationBasis(renderer, parameters, keys, arguments.seconds, seeds)
        gate = ObservationEnergyGate(
            basis, reference, renderer.sample_rate, arguments.tolerance_db
        )
        measurement = TorchMetallicLoss(loss)
        validation = [measurement.validate(a) for a, _ in basis.bases.values()]
        weights, solver = gated_weights(basis, measurement, gate, arguments.iterations)
        values = dict(parameters, **dict(zip(keys, (20 * np.log10(weights)).tolist())))
        actual = [renderer.render(values, arguments.seconds, seed) for seed in seeds]
        errors = [gate.actual_errors(audio).tolist() for audio in actual]
        if (
            not np.isfinite(errors).all()
            or np.max(np.abs(errors)) > arguments.tolerance_db + 0.02
        ):
            raise ValueError("Fresh actual render failed the early-energy gate")
        arguments.output.mkdir(parents=True, exist_ok=True)
        write_wav(
            arguments.output / "reference.wav",
            AudioBuffer(reference, renderer.sample_rate),
        )
        search = Search(
            renderer,
            loss,
            arguments.output,
            arguments.seconds,
            arguments.target.title(),
            seeds,
        )
        search.parameters = values
        search.history = [
            dict(
                stage="full shared loss with per-seed early RMS gates",
                **solver,
                bounds_db={key: [-45, 6] for key in keys},
                tolerance_db=arguments.tolerance_db,
                bins_seconds=gate.bins,
                actual_seed_errors_db=errors,
                basis_validation=basis.validation,
                measurement_validation=validation,
                fixed_parameters={
                    key: value for key, value in parameters.items() if key not in keys
                },
            )
        ]
        search.save()
        verify_candidate(renderer, arguments.output)
        print(json.dumps(search.history[-1]), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("crash", "ride", "gong", "hihat"))
    parser.add_argument("start", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seconds", type=float, default=12)
    parser.add_argument("--iterations", type=int, default=25)
    parser.add_argument("--tolerance-db", type=float, default=3.5)
    refine(parser.parse_args())
