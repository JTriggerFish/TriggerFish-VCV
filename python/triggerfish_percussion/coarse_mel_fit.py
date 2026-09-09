"""Perceptual refinement of broad observation curves with an onset/bloom guard."""

import json
import numpy as np
from scipy.optimize import minimize
import torch

from .coarse_observation_fit import CoarseObservationBasis
from .perceptual_observation_loss import TorchAuralossMel
from .spectral_bloom_basis import SpectralBloomGuard


def polish_coarse_mel(search, mel, knots=(120, 600, 3000, 15000)):
    incoming = dict(search.parameters)
    incoming_audio = [search.audio(incoming, seed) for seed in search.seeds]
    before = float(np.mean([mel.score(audio) for audio in incoming_audio]))
    basis = CoarseObservationBasis(
        search.renderer, search.parameters, search.seconds, search.seeds, knots
    )
    start = np.linalg.lstsq(
        basis.weights,
        10 ** (np.array([search.parameters[k] for k in basis.keys]) / 20),
        rcond=None,
    )[0]
    start = np.clip(start, 10 ** (-45 / 20), 10 ** (6 / 20))
    # Re-centre at the projected coarse curve for optimization. The guard and
    # acceptance comparison separately protect the actual incoming audio.
    basis.bases = {
        seed: (audio + (start - basis.amplitudes) @ columns, columns)
        for seed, (audio, columns) in basis.bases.items()
    }
    basis.amplitudes = start
    guard = SpectralBloomGuard(basis, search.loss, 1.0, baseline_audio=incoming_audio)
    measurement = TorchAuralossMel(mel)
    tensors = [(torch.tensor(a), torch.tensor(c)) for a, c in basis.bases.values()]
    for audio, _ in basis.bases.values():
        measurement.validate(audio)
    origin = torch.tensor(start)
    candidates = []

    def objective(values):
        weights = torch.tensor(values, requires_grad=True)
        loss = sum(measurement(a + (weights - origin) @ c) for a, c in tensors) / len(
            tensors
        )
        loss.backward()
        value = float(loss.detach())
        gradient = weights.grad.numpy().copy()
        if not np.isfinite(value) or not np.isfinite(gradient).all():
            raise ValueError("Nonfinite coarse perceptual objective or gradient")
        candidates.append((value, np.array(values)))
        return value, gradient

    projected_before = objective(start)[0]
    result = minimize(
        objective,
        start,
        jac=True,
        method="SLSQP",
        bounds=[(10 ** (-45 / 20), 10 ** (6 / 20))] * len(knots),
        constraints=[
            dict(type="ineq", fun=guard.values, jac=lambda x: guard.evaluate(x)[1])
        ],
        options=dict(maxiter=40, ftol=1e-8),
    )
    feasible = [row for row in candidates if np.min(guard.values(row[1])) >= -1e-5]
    chosen, actual, selected = start, before, False
    if feasible:
        after, chosen = min(feasible, key=lambda row: row[0])
        candidate = basis.parameters(chosen)
        actual = float(
            np.mean([mel.score(search.audio(candidate, seed)) for seed in search.seeds])
        )
        if abs(actual - after) > 1e-5:
            raise ValueError(
                f"Actual coarse perceptual score differs from basis: {actual} vs {after}"
            )
        selected = actual < before
        if selected:
            search.parameters = candidate
    search.history.append(
        dict(
            stage=f"{len(knots)}-coordinate Mel polish",
            knots_hz=knots,
            knots_db=(20 * np.log10(chosen)).tolist(),
            before=before,
            projected_before=projected_before,
            candidate_score=actual if feasible else None,
            after=float(actual) if selected else before,
            selected=selected,
            guard=guard.specification,
            basis_validation=basis.validation,
            solver_message=str(result.message),
            iterations=int(result.nit),
        )
    )
    search.save()
    print(json.dumps(search.history[-1]), flush=True)
