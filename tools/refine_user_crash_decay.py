"""Preserve a user crash's front while testing existing damping/transport controls.

No DSP changes, new knots, bar edits, output gain matching or publication.
The saved user's gesture is archived separately; fitting uses the reference cell.
"""

import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import least_squares

from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


class ProtectedDecay:
    """Equal-weight early/late decay regions, plus a user-front preservation guard."""

    def __init__(self, reference, original, rate):
        self.decay = BandDecayShapeLoss(reference, rate)
        self.front = SpectralBloomLoss(original, rate)
        self.front_mask = self.front.active
        self.regions = [
            (self.decay.times >= a) & (self.decay.times < b)
            for a, b in ((0.5, 1), (1, 2.5), (2.5, 5.5))
        ]
        self.specification = dict(
            decay=self.decay.specification,
            front_tolerance_db=1.5,
            front_seconds=0.45,
            region_seconds=[[0.5, 1], [1, 2.5], [2.5, 5.5]],
        )

    def residual(self, audio):
        difference = self.decay.relative_db(self.decay.power(audio)) - self.decay.target
        parts = []
        for region in self.regions:
            band_errors = [
                error[mask & region]
                for error, mask in zip(difference, self.decay.mask)
                if np.count_nonzero(mask & region) >= 4
            ]
            parts.extend(
                error / np.sqrt(error.size * len(band_errors)) for error in band_errors
            )
        front = (self.front.db(self.front.power(audio)) - self.front.target)[
            self.front_mask, :5
        ]
        violation = np.maximum(np.abs(front) - 1.5, 0)
        return np.r_[
            np.concatenate(parts) / np.sqrt(3),
            3 * violation.ravel() / np.sqrt(violation.size),
        ]

    def diagnostics(self, audio):
        front = (self.front.db(self.front.power(audio)) - self.front.target)[
            self.front_mask, :5
        ]
        return dict(
            decay_db=self.decay.diagnostics(audio)["shape_error_db"],
            front_rms_db=float(np.sqrt(np.mean(front**2))),
            score=float(np.linalg.norm(self.residual(audio))),
        )


def trial(renderer, base, objective, reference, output, coupled):
    keys = ["body_decay_seconds_0", "body_decay_seconds_7"]
    low, high = [0.1, 0.1], [30, 10]
    if coupled:
        keys += ["bloom_rate", "bloom_energy_acceleration"]
        low += [0.2, 0.001]
        high += [16, 0.5]
    origin = np.log([base[k] for k in keys])
    if coupled:
        origin[-1] = np.log(0.1)
    records, cache = [], {}

    def evaluate(x):
        key = tuple(x)
        if key not in cache:
            parameters = dict(
                base, **dict(zip(keys, np.clip(np.exp(x), low, high).tolist()))
            )
            audio = renderer.render(parameters, 6)
            residual = objective.residual(audio)
            cache[key] = residual
            records.append(
                dict(parameters=parameters, score=float(np.linalg.norm(residual)))
            )
            if len(records) % 25 == 0:
                print(
                    json.dumps(
                        dict(
                            coupled=coupled,
                            evaluations=len(records),
                            best=min(r["score"] for r in records),
                        )
                    ),
                    flush=True,
                )
        return cache[key]

    def jacobian(x):
        columns = []
        for i in range(len(x)):
            a, b = x.copy(), x.copy()
            a[i] = max(np.log(low[i]), x[i] - 0.02)
            b[i] = min(np.log(high[i]), x[i] + 0.02)
            columns.append((evaluate(b) - evaluate(a)) / (b[i] - a[i]))
        return np.array(columns).T

    evaluate(np.log([base[k] for k in keys]))
    result = least_squares(
        evaluate,
        origin,
        jac=jacobian,
        bounds=(np.log(low), np.log(high)),
        max_nfev=28,
        xtol=0.003,
        ftol=0.003,
        gtol=0.001,
        x_scale="jac",
    )
    chosen = min(records, key=lambda r: r["score"])
    audio = renderer.render(chosen["parameters"], 6)
    checkpoint(
        renderer,
        objective.front,
        output,
        "Candidate crash — decay refinement",
        chosen["parameters"],
        reference,
        [
            dict(
                stage="protected-front decay",
                specification=objective.specification,
                controls=keys,
                bounds=list(zip(low, high)),
                solver=str(result.message),
                log_probe=0.02,
                trace=records,
                diagnostics=objective.diagnostics(audio),
            )
        ],
    )
    print(
        json.dumps(
            dict(
                output=str(output),
                **objective.diagnostics(audio),
                controls={k: chosen["parameters"][k] for k in keys},
            )
        ),
        flush=True,
    )


def main():
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    output = Path("build/crash-user-decay")
    try:
        source = verify_candidate(renderer, output / "before")
        base = source["parameters"]
        reference = aligned_reference(renderer, 6)
        original = renderer.render(base, 6)
        objective = ProtectedDecay(reference, original, renderer.sample_rate)
        print(json.dumps(dict(baseline=objective.diagnostics(original))), flush=True)
        for coupled in (False, True):
            trial(
                renderer,
                base,
                objective,
                reference,
                output / ("coupled" if coupled else "damping-only"),
                coupled,
            )
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
