"""Restore the approved low pair without undoing the upper-texture refinement.

Only the 120/240 Hz observation amplitudes change. Match their early output
energy to the approved preset, accounting for changed excitation/diffusion.
This is a sound-preservation correction, not another reference-target fit.
"""

import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import least_squares
from scipy.signal import stft

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.observation_fit_basis import ObservationBasis
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def low_levels(audio, rate):
    f, t, z = stft(audio, rate, nperseg=8192, noverlap=8192 - round(0.01 * rate))
    power = abs(z) ** 2
    return np.array(
        [
            10
            * np.log10(
                max(
                    1e-15,
                    power[(f >= lo) & (f < hi)][:, (t >= 0.05) & (t < 0.5)]
                    .sum(axis=0)
                    .mean(),
                )
            )
            for lo, hi in ((80, 180), (180, 300))
        ]
    )


def run():
    torch.set_num_threads(1)
    root = Path.cwd()
    output = root / "build/gong-low-body-restored"
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", root)
    try:
        approved = json.loads(
            (root / "build/gong-layered-timing/baseline/search.json").read_text()
        )["parameters"]
        current = dict(r.initial)
        keys = ["resolved_level_0", "resolved_level_1"]
        if [current[f"resolved_frequency_{i}"] for i in range(2)] != [120, 240]:
            raise ValueError("Expected the approved 120/240 Hz pair")
        base = dict(current, **{key: approved[key] for key in keys})
        seeds = (1675, 1982)
        targets = np.array(
            [low_levels(r.render(approved, 6, seed), r.sample_rate) for seed in seeds]
        )
        basis = ObservationBasis(r, base, keys, 6, seeds)

        def residual(levels):
            p = dict(base, **dict(zip(keys, levels)))
            db = np.array(
                [low_levels(basis.render(p, 6, seed), r.sample_rate) for seed in seeds]
            )
            return (db - targets).ravel()

        result = least_squares(
            residual,
            [base[k] for k in keys],
            bounds=([-30, -30], [0, 0]),
            diff_step=0.01,
        )
        p = dict(base, **{k: float(v) for k, v in zip(keys, result.x)})
        ref = aligned_reference(r, 6)
        loss = LayeredBandLoss(ref, r.sample_rate, audibility=True)
        checkpoint(
            r,
            loss,
            output,
            "Gong — restored low body and metallic sizzle",
            p,
            ref,
            [
                dict(
                    stage="preserve approved low body",
                    keys=keys,
                    bands_hz=[[80, 180], [180, 300]],
                    region_seconds=[0.05, 0.5],
                    approved_target_db=targets.tolist(),
                    errors_db=residual(result.x).tolist(),
                    basis_validation=basis.validation,
                )
            ],
        )
        verify_candidate(r, output)
        report = dict(
            levels={k: [current[k], p[k]] for k in keys},
            target_db=targets.tolist(),
            actual_db=[
                low_levels(r.render(p, 6, s), r.sample_rate).tolist() for s in seeds
            ],
            changed_keys=[k for k in p if p[k] != current[k]],
        )
        (output / "restoration.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    run()
