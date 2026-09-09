"""Checkpoint a texture-only gong candidate for render/provenance review."""

import json
import os
from pathlib import Path

import numpy as np
import torch
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def main():
    torch.set_num_threads(1)
    output = Path("build/gong-beating-recovery")
    screen = json.loads((output / "screen.json").read_text(encoding="utf8"))
    chosen = next(r for r in screen["variants"] if r["name"] == "clean-low-0.4")
    base = next(r["parameters"] for r in screen["variants"] if r["name"] == "before")
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        ref = aligned_reference(r, 6)
        loss = SpectralBloomLoss(ref, r.sample_rate)
        old = r.render(base, 6)
        new = r.render(chosen["parameters"], 6)
        front = SpectralBloomLoss(old, r.sample_rate)
        delta = (front.db(front.power(new)) - front.target)[front.active, :5]
        history = [
            dict(
                stage="texture recovery",
                screen="screen.json",
                selected=chosen["name"],
                front_change_db=float(np.sqrt(np.mean(delta**2))),
                fixed="modal series, damping, bloom, observation, strike and gain",
            )
        ]
        checkpoint(
            r, loss, output / "before", "Gong — before texture recovery", base, ref, []
        )
        checkpoint(
            r,
            loss,
            output / "candidate",
            "Gong — gentler beating",
            chosen["parameters"],
            ref,
            history,
        )
        verify_candidate(r, output / "candidate")
        print(json.dumps(history))
    finally:
        r.close()


if __name__ == "__main__":
    main()
