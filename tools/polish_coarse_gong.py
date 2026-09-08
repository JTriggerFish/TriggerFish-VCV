"""Library Mel refinement of four broad observation coordinates, with bloom guard."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from triggerfish_percussion.coarse_observation_fit import CoarseObservationBasis
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.perceptual_observation_loss import TorchAuralossMel
from triggerfish_percussion.spectral_bloom_basis import SpectralBloomGuard
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def polish(search, mel):
    knots = (120, 600, 3000, 15000)
    basis = CoarseObservationBasis(
        search.renderer, search.parameters, 6, search.seeds, knots
    )
    start = np.linalg.lstsq(
        basis.weights,
        10 ** (np.array([search.parameters[k] for k in basis.keys]) / 20),
        rcond=None,
    )[0]
    start = np.clip(start, 10 ** (-45 / 20), 10 ** (6 / 20))
    # Re-centre the same exact affine basis at the actual incoming curve so
    # the bloom guard protects that curve, not the uniform probing baseline.
    basis.bases = {
        seed: (audio + (start - basis.amplitudes) @ columns, columns)
        for seed, (audio, columns) in basis.bases.items()
    }
    basis.amplitudes = start
    guard = SpectralBloomGuard(basis, search.loss, 1.0)
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

    before = objective(start)[0]
    result = minimize(
        objective,
        start,
        jac=True,
        method="SLSQP",
        bounds=[(10 ** (-45 / 20), 10 ** (6 / 20))] * 4,
        constraints=[
            dict(type="ineq", fun=guard.values, jac=lambda x: guard.evaluate(x)[1])
        ],
        options=dict(maxiter=40, ftol=1e-8),
    )
    feasible = [row for row in candidates if np.min(guard.values(row[1])) >= -1e-5]
    after, chosen = min(feasible, key=lambda row: row[0])
    search.parameters = basis.parameters(chosen)
    actual = np.mean(
        [mel.score(search.audio(search.parameters, seed)) for seed in search.seeds]
    )
    if abs(actual - after) > 1e-5:
        raise ValueError("Actual coarse perceptual score differs from basis")
    search.history.append(
        dict(
            stage="four-coordinate Mel polish",
            knots_hz=knots,
            knots_db=(20 * np.log10(chosen)).tolist(),
            before=before,
            after=float(actual),
            guard=guard.specification,
            basis_validation=basis.validation,
            solver_message=str(result.message),
            iterations=int(result.nit),
        )
    )
    search.save()
    print(json.dumps(search.history[-1]), flush=True)


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
        search = checkpoint(
            renderer,
            shape,
            args.output,
            "Gong protected harmonic core",
            saved["parameters"],
            reference,
            saved["history"] + [dict(parent=str(args.source))],
        )
        search.seeds = (1675, 1776)
        polish(search, AuralossMel(reference, renderer.sample_rate))
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
