"""Compare losses in the same exact-render observation subspace; publish nothing.

This is a controlled optimizer experiment, not a complete instrument fitter.
The fine-spectrum residual is experimental, not a perceptual acceptance score.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.regional_spectrum_audit import RegionalSpectrumAudit
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def run(args):
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], args.target, Path.cwd())
    try:
        saved = verify_candidate(r, args.source)
        p = saved["parameters"]
        if len(set(args.keys)) != len(args.keys) or any(
            key not in p or not -60 <= p[key] <= 6 for key in args.keys
        ):
            raise ValueError("Expected unique active observation keys within -60..6 dB")
        seed = int(r.metadata["event"]["seed"])
        reference = read_wav(args.source / "reference.wav").mono().samples
        bands = LayeredBandLoss(reference, r.sample_rate, audibility=True)
        spectrum = RegionalSpectrumAudit(
            reference,
            r.sample_rate,
            ((0.04, 0.2), (0.2, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 6)),
        )
        basis = ObservationBasis(r, p, args.keys, 6, (seed,))

        def parameters(db):
            return dict(p, **dict(zip(args.keys, map(float, db))))

        def residual(db, detailed):
            audio = basis.render(parameters(db), 6, seed)
            if detailed:
                error = spectrum.db(spectrum.power(audio)) - spectrum.target
                # Fixed reference mask; no protected frequencies or bar levels.
                return error[spectrum.active]
            return bands.residual(bands.envelopes(audio)).ravel()

        rows = []
        for detailed in (False, True):
            for start in (
                np.array([p[k] for k in args.keys]),
                np.full(len(args.keys), -24.0),
            ):
                result = least_squares(
                    lambda db: residual(db, detailed),
                    start,
                    bounds=(-60, 6),
                    diff_step=0.001,
                    max_nfev=80,
                    ftol=1e-8,
                    xtol=1e-6,
                )
                candidate = parameters(result.x)
                actual = r.render(candidate, 6, seed)
                predicted = basis.render(candidate, 6, seed)
                maximum = float(np.max(abs(actual - predicted)))
                if not np.isfinite(maximum) or maximum > 3e-5:
                    raise ValueError("Candidate failed exact-render basis validation")
                row = dict(
                    objective=(
                        "fine spectrum (experimental)" if detailed else "layered bands"
                    ),
                    start_db=start.tolist(),
                    levels_db=result.x.tolist(),
                    parameters=candidate,
                    solver_success=bool(result.success),
                    evaluations=result.nfev,
                    exact_render_max_error=maximum,
                    layered=bands.attribution(actual),
                    spectrum=spectrum.measure(actual),
                )
                rows.append(row)
                print(
                    json.dumps(
                        {k: row[k] for k in ("objective", "levels_db", "evaluations")}
                    ),
                    flush=True,
                )
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "objective-profile.json").write_text(
            json.dumps(
                dict(
                    source=str(args.source),
                    metadata=r.metadata,
                    keys=args.keys,
                    basis_validation=basis.validation,
                    rows=rows,
                ),
                indent=2,
            )
        )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--keys", nargs="+", required=True)
    run(parser.parse_args())
