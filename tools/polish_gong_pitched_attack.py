"""Restore whole-gong bloom around a cleaner, explicitly measured low core.

The first-half-second low-band MR-STFT is the objective. Hinge penalties guard
whole-sound Mel (+5%), band envelopes (+10%) and bloom rise (+10%) relative to
the ORIGINAL preset, not the warm start. Recheck all metrics across held-out
seeds before publishing. No mode frequencies/levels, damping or gains are fitted.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.signal import butter, sosfilt

from triggerfish_percussion.attack_ridge_loss import AttackRidgeLoss
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


class PitchedAttackLoss:
    units = "normalized low attack MR-STFT plus whole-sound regression penalties"

    def __init__(self, reference, rate, baselines):
        self.lowpass = butter(6, 1500, fs=rate, output="sos")
        self.attack = AttackRidgeLoss(sosfilt(self.lowpass, reference), rate, 0.5)
        self.mel = ReferenceFloorMel(reference, rate, 60)
        self.bloom = SpectralBloomLoss(reference, rate)
        self.baseline = np.maximum(
            np.mean([self.measure(x) for x in baselines], axis=0), 1e-12
        )
        self.limits = self.baseline[1:] * np.array([1.05, 1.10, 1.10])
        self.specification = dict(
            attack=self.attack.specification,
            lowpass_hz=1500,
            baseline=self.baseline.tolist(),
            limits=self.limits.tolist(),
            hinge_weight=10,
            normalization=False,
            metric_order=["attack", "mel", "envelope", "rise"],
        )

    def measure(self, audio):
        b = self.bloom.diagnostics(audio)
        return np.array(
            [
                self.attack.score(sosfilt(self.lowpass, audio)),
                self.mel.score(audio),
                b["envelope_rms_db"],
                b["rise_rms_db"],
            ]
        )

    def score(self, audio):
        m = self.measure(audio)
        return float(
            m[0] / self.baseline[0] + 10 * np.maximum(m[1:] / self.limits - 1, 0).sum()
        )


def main(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        rows = json.loads((args.screen / "screen.json").read_text(encoding="utf8"))
        row = next(x for x in rows if x["name"] == args.candidate)
        seeds = tuple(r.metadata["event"]["seed"] + i for i in (0, 307))
        ref = aligned_reference(r, 6)
        baseline = [r.render(rows[0]["parameters"], 6, s) for s in seeds]
        loss = PitchedAttackLoss(ref, r.sample_rate, baseline)
        search = checkpoint(
            r,
            loss,
            args.output,
            "Gong — clearer pitched attack",
            row["parameters"],
            ref,
            [
                dict(
                    stage="fixed five-mode low core",
                    source=str(args.screen),
                    trial=args.candidate,
                )
            ],
        )
        search.seeds = seeds
        refine_scalar(
            search,
            dict(
                bloom_rate=(0.5, 6),
                bloom_energy_acceleration=(0, 0.3),
                bloom_energy_sensitivity=(0, 0.6),
                body_brightness=(-65, -25),
            ),
            budget=args.budget,
            step=0.01,
            method="Powell",
        )
        verify_candidate(r, args.output)
        audit = []
        for offset in (0, 307, 911, 1601):
            seed = r.metadata["event"]["seed"] + offset
            audit.append(
                dict(
                    seed=seed,
                    before=loss.measure(
                        r.render(rows[0]["parameters"], 6, seed)
                    ).tolist(),
                    after=loss.measure(r.render(search.parameters, 6, seed)).tolist(),
                )
            )
        (args.output / "audit.json").write_text(
            json.dumps(audit, indent=2), encoding="utf8"
        )
        print(json.dumps(dict(audit=audit)), flush=True)
    finally:
        r.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("screen", type=Path)
    p.add_argument("candidate")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--budget", type=int, default=120)
    main(p.parse_args())
