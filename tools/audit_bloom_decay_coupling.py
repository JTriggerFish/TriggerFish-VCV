"""Measure how far two visible damping endpoints can compensate a bloom edit.

This is an offline control-surface experiment, not output normalization or a
new runtime voice. All steps use the current exact Wasm and fixed gesture.
"""

import json
import os
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss


def main():
    output = Path("build/crash-user-decay")
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        p = verify_candidate(r, output / "before")["parameters"]
        loss = SpectralBloomLoss(r.render(p, 6), r.sample_rate)
        # Preserve the incoming patch's audible late envelope; don't chase its
        # silence floor, and don't normalize either input/output audio.
        mask = np.zeros_like(loss.target, dtype=bool)
        mask[:, 9:] = loss.target[:, 9:] > loss.target.max() - 55

        def measurement(parameters):
            return (loss.db(loss.power(r.render(parameters, 6))) - loss.target)[mask]

        keys = ["body_decay_seconds_0", "body_decay_seconds_7"]
        columns = []
        for key in keys:
            a, b = dict(p), dict(p)
            a[key] = max(0.1, p[key] / 1.04)
            b[key] = min(30, p[key] * 1.04)
            columns.append((measurement(b) - measurement(a)) / np.log(b[key] / a[key]))
        jacobian = np.array(columns).T
        changed = dict(p, bloom_rate=p["bloom_rate"] * 1.25)
        before = measurement(changed)
        result = lsq_linear(
            jacobian,
            -before,
            bounds=(
                np.log(0.1 / np.array([p[k] for k in keys])),
                np.log(30 / np.array([p[k] for k in keys])),
            ),
        )
        compensated = dict(
            changed,
            **{
                k: float(np.clip(p[k] * np.exp(x), 0.1, 30))
                for k, x in zip(keys, result.x)
            },
        )
        after = measurement(compensated)
        report = dict(
            edit="diffusion rate +25%",
            region="1–6 s; above incoming patch peak −55 dB",
            probe="±4% in damping time; log coordinates",
            jacobian_singular_values=np.linalg.svd(jacobian, compute_uv=False).tolist(),
            uncompensated_rms_db=float(np.sqrt(np.mean(before**2))),
            compensated_rms_db=float(np.sqrt(np.mean(after**2))),
            predicted_rms_db=float(
                np.sqrt(np.mean((before + jacobian @ result.x) ** 2))
            ),
            compensating_controls={k: compensated[k] for k in keys},
            parameter_change={k: compensated[k] / p[k] for k in keys},
            damping_only_limits="nonnegative loss, T60 <=30s; no gain restoration",
        )
        (output / "compensation-audit.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report), flush=True)
    finally:
        r.close()


if __name__ == "__main__":
    main()
