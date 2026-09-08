"""Fit metallic observation amplitudes to absolute band/time energy.

Development-only diagnostic subproblem; never publishes a workbench preset.
Run from the repository with the analysis environment and EMSDK_NODE configured.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.regional_energy_basis import RegionalEnergyBasis
from triggerfish_percussion.regional_energy_loss import RegionalEnergyLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def fit(renderer, source, output, specification, iterations):
    """Hold all non-observation controls and the recorded strike fixed."""
    saved = verify_candidate(renderer, source)
    seconds = saved["duration_seconds"]
    reference = aligned_reference(renderer, seconds)
    loss = RegionalEnergyLoss(reference, renderer.sample_rate, **specification)
    search = Search(renderer, loss, output, seconds, "Regional observation trial")
    search.parameters = saved["parameters"]
    polish(search, iterations)
    write_wav(output / "reference.wav", AudioBuffer(reference, renderer.sample_rate))
    search.save()
    verify_candidate(renderer, output)
    print(json.dumps(loss.diagnostics(search.audio(search.parameters))))


def polish(search, iterations):
    """Reusable exact-energy observation solve, including every training seed."""
    loss, renderer, seconds = search.loss, search.renderer, search.seconds
    keys = [
        key
        for key, value in search.parameters.items()
        if key.startswith("resolved_level_") and value > -71.99
    ]
    basis = ObservationBasis(renderer, search.parameters, keys, seconds, search.seeds)
    energy = RegionalEnergyBasis(basis, loss)
    # Use the same positive-amplitude bounds as the full-spectrum fitter.
    low, high = 10 ** (-45 / 20), 10 ** (6 / 20)
    result = least_squares(
        lambda x: energy.evaluate(x)[0],
        np.clip(basis.amplitudes, low, high),
        jac=lambda x: energy.evaluate(x)[1],
        bounds=(low, high),
        max_nfev=iterations,
    )
    validation = energy.validate(basis, result.x)
    candidate = dict(
        search.parameters, **dict(zip(keys, (20 * np.log10(result.x)).tolist()))
    )
    before = float(np.linalg.norm(search.residual(search.parameters)))
    after = float(np.linalg.norm(search.residual(candidate)))
    if after < before:
        search.parameters = candidate
    search.history.append(
        dict(
            stage="exact quadratic regional observation",
            before=before,
            after=after,
            selected=after < before,
            basis_validation=basis.validation,
            quadratic_error=validation,
            solver_message=result.message,
            evaluations=int(result.nfev),
            bounds_db=[-45, 6],
            parameters=dict(search.parameters),
        )
    )
    search.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("crash", "ride", "gong", "hihat"))
    parser.add_argument("source", type=Path, help="Verified candidate directory")
    parser.add_argument("output", type=Path, help="New diagnostic directory")
    parser.add_argument(
        "measurement", type=Path, help="JSON with bands and regions, in Hz and seconds"
    )
    parser.add_argument("--iterations", type=int, default=300)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("Iterations must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        fit(
            renderer,
            args.source,
            args.output,
            json.loads(args.measurement.read_text(encoding="utf8")),
            args.iterations,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
