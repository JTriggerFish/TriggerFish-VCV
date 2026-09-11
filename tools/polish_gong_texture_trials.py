"""Polish a listening family using only shared series gain and tilt.

No EQ, frequency movement, extra decay knots, hard gain curve, or bass lock.
Texture-family differences remain explicit and are not optimized away.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import LinearConstraint, minimize

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.regional_observation_basis import RegionalObservationBasis
from triggerfish_percussion.structured_modal_levels import StructuredModalLevels
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from gong_texture_comparison import GongTextureComparison


def polish(r, target, p, output):
    surface = StructuredModalLevels(p)
    keys = surface.keys
    print(
        json.dumps(dict(stage="basis", output=str(output), handles=len(keys))),
        flush=True,
    )
    basis = ObservationBasis(r, p, keys, 6, (1675,))
    cache = RegionalObservationBasis(basis, target.spectrum, 1675)
    cache.validate(basis, target.spectrum, 1675, basis.amplitudes * 0.8)

    def objective(controls):
        db = surface.levels(controls)
        power = cache.power(10 ** (db / 20))
        spectral = target.power_residual(power)
        return float(
            np.mean(spectral**2) + 0.15**2 * np.mean((db - surface.initial) ** 2)
        )

    result = minimize(
        objective,
        np.zeros(2),
        method="SLSQP",
        bounds=((-12, 12), (-3, 3)),
        constraints=LinearConstraint(
            surface.matrix,
            -71.989 - surface.initial,
            6 - surface.initial,
        ),
        options=dict(maxiter=150, ftol=1e-8, eps=1e-3),
    )
    levels = surface.levels(result.x)
    if not result.success or np.any(levels <= -71.99) or np.any(levels > 6 + 1e-7):
        raise ValueError(f"Shared series fit failed: {result.message}")
    fitted = dict(p, **dict(zip(keys, map(float, levels))))
    cache_error = cache.validate(basis, target.spectrum, 1675, 10 ** (levels / 20))
    actual = r.render(fitted, 6, 1675)
    difference = float(np.max(abs(actual - basis.render(fitted, 6, 1675))))
    if not np.isfinite(difference) or difference > 3e-5:
        raise ValueError("Exact candidate differs from the observation basis")
    ref = aligned_reference(r, 6)
    measurements = [target.measure(r.render(fitted, 6, seed)) for seed in (1675, 1982)]
    history = [
        dict(
            stage="EQ-free texture comparison",
            target=target.specification,
            regularization=".15 times RMS shared series level change",
            fit_coordinates=["series_gain_db", "series_tilt_db_per_octave"],
            coordinate_values=list(map(float, result.x)),
            affected_keys=keys,
            evaluations=result.nfev,
            solver_success=bool(result.success),
            cache_error_db=cache_error,
            exact_render_max_error=difference,
            measures=measurements,
        )
    ]
    checkpoint(r, target, output, "Gong texture trial", fitted, ref, history)
    verify_candidate(r, output)
    print(json.dumps(dict(output=str(output), measurements=measurements)), flush=True)
    return fitted


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = json.loads((args.directory / "baseline.json").read_text())["parameters"]
        ref = aligned_reference(r, 6)
        loss = GongTextureComparison(ref, r.render(base, 6, 1675), r.sample_rate)
        rows = json.loads((args.directory / "screen.json").read_text())
        for family in args.families:
            choices = [x for x in rows if x["family"] == family]
            winner = min(choices, key=lambda x: x["score"])
            polish(r, loss, winner["parameters"], args.directory / family)
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--families", nargs="+", default=["movement"])
    run(parser.parse_args())
