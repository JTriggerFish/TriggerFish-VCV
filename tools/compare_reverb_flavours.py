#!/usr/bin/env python3
"""Compare Base and Optimized coefficients on the current architecture."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from triggerfish_reverb.objectives import aggregate_resonance, resonance_per_control
from triggerfish_reverb.velvet import DifferentiableVelvetReverb

from optimize_velvet_reverb import ARTIFACT_FORMAT, control_grid


def load_artifact(model: DifferentiableVelvetReverb, path: Path) -> None:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("format") != ARTIFACT_FORMAT:
        raise ValueError("artifact does not describe the current two-stage VFM")
    reference = model.raw_main_ratios
    model.set_coefficients(
        reference.new_tensor(artifact["main_delay_ratio"]),
        reference.new_tensor(artifact["velvet_delay_ms"]),
        torch.tensor(
            artifact["permutations"],
            dtype=model.permutations.dtype,
            device=reference.device,
        ),
        reference.new_tensor(artifact["signs"]),
    )


@torch.no_grad()
def score(
    model: DifferentiableVelvetReverb,
    frequencies: torch.Tensor,
) -> tuple[float, float]:
    values = []
    worst_band = 0.0
    controls = control_grid()
    for first in range(0, len(controls), 5):
        wall, resolvent = model.response(frequencies, controls[first : first + 5])
        combined = torch.cat((wall.flatten(2), resolvent.flatten(2)), dim=2)
        per_control, diagnostics = resonance_per_control(combined, frequencies)
        values.append(per_control)
        worst_band = max(worst_band, float(diagnostics["worst_band"]))
    return float(aggregate_resonance(torch.cat(values))), worst_band


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--seconds", type=float, default=8.0)
    parser.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu"
    )
    args = parser.parse_args()
    if args.seconds <= 0.0:
        parser.error("seconds must be positive")

    device = torch.device(args.device)
    frequencies = torch.arange(20.0, 20_000.0, 1.0 / args.seconds, device=device)
    baseline = DifferentiableVelvetReverb().to(device=device)
    optimized = DifferentiableVelvetReverb().to(device=device)
    load_artifact(optimized, args.artifact)
    for name, model in (("Base", baseline), ("Optimized", optimized)):
        objective, worst_band = score(model, frequencies)
        print(f"{name}: objective={objective:.9g} worst-band={worst_band:.9g}")


if __name__ == "__main__":
    main()
