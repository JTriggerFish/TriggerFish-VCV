"""Fit exact-render observation amplitudes with optional analysis autograd.

Torch differentiates the analysis and a validated affine combination of actual
C++ renders. It does not differentiate or approximate the nonlinear instrument.
"""

import json
import numpy as np
from scipy.optimize import minimize
import torch

from .metallic_balance_loss import MetallicBalanceLoss
from .observation_fit_basis import ObservationBasis
from .torch_metallic_loss import TorchMetallicLoss
from .perceptual_fit_losses import AuralossMel


def _comparison_score(search, parameters):
    score = float(np.linalg.norm(search.residual(parameters)))
    # ScalarAudioLoss exposes sqrt(score) only as a residual adapter. Reports
    # for the library objective must retain its native units, not their root.
    return score * score if isinstance(search.loss, AuralossMel) else score


def _observation_keys(parameters, fixed_keys):
    if not set(fixed_keys).issubset(parameters) or any(
        not key.startswith("resolved_level_") for key in fixed_keys
    ):
        raise ValueError("Only observation-bar keys can be held fixed here")
    return [
        key
        for key, value in parameters.items()
        if key.startswith("resolved_level_")
        and value > -71.99
        and key not in fixed_keys
    ]


def polish_observation_autograd(
    search, iterations=40, constraint_factory=None, bounds_db=(-45, 6), fixed_keys=()
):
    minimum_db, maximum_db = bounds_db
    if not np.isfinite(bounds_db).all() or not -71.99 < minimum_db < maximum_db <= 6:
        raise ValueError(
            "Observation bounds must retain positive active modes inside the UI range"
        )
    if isinstance(search.loss, MetallicBalanceLoss):
        measurement = TorchMetallicLoss(search.loss)
    elif isinstance(search.loss, AuralossMel):
        from .perceptual_observation_loss import TorchAuralossMel

        measurement = TorchAuralossMel(search.loss)
    else:
        raise ValueError(
            "Analysis autograd requires the metallic balance or auraloss mel loss"
        )
    torch.set_num_threads(1)
    keys = _observation_keys(search.parameters, fixed_keys)
    if not keys:
        return
    before = _comparison_score(search, search.parameters)
    basis = ObservationBasis(
        search.renderer, search.parameters, keys, search.seconds, search.seeds
    )
    guard = constraint_factory(basis) if constraint_factory is not None else None
    validation = [
        measurement.validate(basis.render(search.parameters, search.seconds, seed))
        for seed in search.seeds
    ]
    prepared = {
        seed: (torch.tensor(audio), torch.tensor(columns))
        for seed, (audio, columns) in basis.bases.items()
    }
    initial = torch.tensor(basis.amplitudes)
    evaluated = []

    def objective(amplitudes):
        weights = torch.tensor(amplitudes, requires_grad=True)
        total = 0.0
        for audio, columns in prepared.values():
            value = measurement(audio + (weights - initial) @ columns) / len(prepared)
            value.backward()
            total += float(value.detach())
        gradient = weights.grad.numpy().copy()
        if not np.isfinite(total) or not np.isfinite(gradient).all():
            raise ValueError("Nonfinite observation objective or gradient")
        if guard is not None:
            evaluated.append((total, np.asarray(amplitudes).copy()))
        return total, gradient

    low, high = 10 ** (minimum_db / 20), 10 ** (maximum_db / 20)
    start = np.clip(basis.amplitudes, low, high)
    # Check the actual objective gradient, not only feature values.
    direction = np.random.default_rng(731).normal(size=len(keys))
    direction /= np.linalg.norm(direction)
    probe = np.clip(start, low + 0.002, high - 0.002)
    step = 1e-5
    numeric = (
        objective(probe + step * direction)[0] - objective(probe - step * direction)[0]
    ) / (2 * step)
    analytic = float(objective(probe)[1] @ direction)
    relative = abs(numeric - analytic) / max(abs(numeric), abs(analytic), 1.0)
    if relative > 0.002:
        raise ValueError(f"Analysis gradient check failed: {numeric=}, {analytic=}")
    validation.append(
        dict(gradient_relative_error=relative, numeric=numeric, analytic=analytic)
    )
    iteration = 0

    def progress(values):
        nonlocal iteration
        iteration += 1
        record = dict(
            stage="autograd observation",
            iteration=iteration,
            parameters=dict(
                search.parameters, **dict(zip(keys, (20 * np.log10(values)).tolist()))
            ),
            metadata=search.renderer.metadata,
            duration_seconds=search.seconds,
            training_seeds=search.seeds,
            status="iteration-checkpoint-not-reviewed",
        )
        search.output.mkdir(parents=True, exist_ok=True)
        temporary = search.output / "autograd-progress.pending.json"
        temporary.write_text(json.dumps(record, indent=2), encoding="utf8")
        temporary.replace(search.output / "autograd-progress.json")
        print(json.dumps(dict(stage=record["stage"], iteration=iteration)), flush=True)

    constraints = (
        []
        if guard is None
        else [
            dict(
                type="ineq",
                fun=lambda x: guard.evaluate(x)[0],
                jac=lambda x: guard.evaluate(x)[1],
            )
        ]
    )
    options = dict(maxiter=iterations, ftol=1e-8)
    if guard is None:
        options.update(gtol=1e-5, maxls=20)
    evaluated.clear()  # Validation probes are not solver iterates.
    result = minimize(
        objective,
        start,
        jac=True,
        method="L-BFGS-B" if guard is None else "SLSQP",
        bounds=[(low, high)] * len(keys),
        callback=progress,
        constraints=constraints,
        options=options,
    )
    chosen = result.x
    # SLSQP may finish on an infeasible line-search step. Retain the best actual
    # feasible evaluation instead of discarding earlier valid improvements.
    if guard is not None:
        for _, point in sorted(evaluated, key=lambda row: row[0]):
            if np.min(guard.values(point)) >= -1e-5:
                chosen = point
                break
    candidate = dict(
        search.parameters, **dict(zip(keys, (20 * np.log10(chosen)).tolist()))
    )
    after = _comparison_score(search, candidate)
    feasible = guard is None or bool(np.min(guard.values(chosen)) >= -1e-5)
    selected = after < before and feasible
    if selected:
        search.parameters = candidate
    search.history.append(
        dict(
            stage="exact-render autograd observation polish",
            before=before,
            after=after,
            score_units=getattr(search.loss, "units", "dB"),
            selected=selected,
            constraint_feasible=feasible,
            constraints=None if guard is None else guard.specification,
            measurement_validation=validation,
            basis_validation=basis.validation,
            iterations=int(result.nit),
            evaluations=int(result.nfev),
            solver_message=str(result.message),
            bounds_db={key: (minimum_db, maximum_db) for key in keys},
            coordinate_scale="linear amplitude",
            fixed_parameters={
                key: value for key, value in basis.initial.items() if key not in keys
            },
            parameters=dict(search.parameters),
        )
    )
    search.save()
    return search.history[-1]
