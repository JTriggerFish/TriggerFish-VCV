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


def polish_observation_autograd(search, iterations=40):
    if not isinstance(search.loss, MetallicBalanceLoss):
        raise ValueError(
            "Analysis autograd requires the validated metallic balance loss"
        )
    torch.set_num_threads(1)
    keys = [
        key
        for key, value in search.parameters.items()
        if key.startswith("resolved_level_") and value > -71.99
    ]
    if not keys:
        return
    before = float(np.linalg.norm(search.residual(search.parameters)))
    basis = ObservationBasis(
        search.renderer, search.parameters, keys, search.seconds, search.seeds
    )
    measurement = TorchMetallicLoss(search.loss)
    validation = [
        measurement.validate(basis.render(search.parameters, search.seconds, seed))
        for seed in search.seeds
    ]
    prepared = {
        seed: (torch.tensor(audio), torch.tensor(columns))
        for seed, (audio, columns) in basis.bases.items()
    }
    initial = torch.tensor(basis.amplitudes)

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
        return total, gradient

    low, high = 10 ** (-45 / 20), 10 ** (6 / 20)
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

    result = minimize(
        objective,
        start,
        jac=True,
        method="L-BFGS-B",
        bounds=[(low, high)] * len(keys),
        callback=progress,
        options=dict(maxiter=iterations, ftol=1e-8, gtol=1e-5, maxls=20),
    )
    candidate = dict(
        search.parameters, **dict(zip(keys, (20 * np.log10(result.x)).tolist()))
    )
    after = float(np.linalg.norm(search.residual(candidate)))
    selected = after < before
    if selected:
        search.parameters = candidate
    search.history.append(
        dict(
            stage="exact-render autograd observation polish",
            before=before,
            after=after,
            selected=selected,
            measurement_validation=validation,
            basis_validation=basis.validation,
            iterations=int(result.nit),
            evaluations=int(result.nfev),
            solver_message=str(result.message),
            bounds_db={key: (-45, 6) for key in keys},
            coordinate_scale="linear amplitude",
            fixed_parameters={
                key: value for key, value in basis.initial.items() if key not in keys
            },
            parameters=dict(search.parameters),
        )
    )
    search.save()
    return search.history[-1]
