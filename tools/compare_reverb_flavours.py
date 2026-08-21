#!/usr/bin/env python3
"""Compare Base and Optimized coefficients on the current architecture."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from triggerfish_reverb.objectives import aggregate_resonance, resonance_per_control
from triggerfish_reverb.velvet import (
    DifferentiableVelvetReverb,
    VelvetControls,
)

from optimize_velvet_reverb import ARTIFACT_FORMAT, control_grid


def validation_control_grid() -> tuple[VelvetControls, ...]:
    """Intermediate controls excluded from the optimizer's training grid."""
    geometries = (
        (0.20, 0.20),
        (0.20, 0.80),
        (0.65, 0.50),
        (0.80, 0.15),
        (0.80, 0.85),
    )
    return tuple(
        VelvetControls(
            space=space,
            aspect=aspect,
            decay=decay,
            damping=damping,
            diffusion=diffusion,
        )
        for space, aspect in geometries
        for decay in (0.25, 0.78)
        for damping in (0.08, 0.58)
        for diffusion in (0.35, 0.90)
    )


def secondary_validation_control_grid() -> tuple[VelvetControls, ...]:
    """A second disjoint validation grid used when selecting a candidate."""
    geometries = (
        (0.35, 0.30),
        (0.35, 0.70),
        (0.90, 0.50),
        (0.90, 0.25),
        (0.90, 0.75),
    )
    return tuple(
        VelvetControls(
            space=space,
            aspect=aspect,
            decay=decay,
            damping=damping,
            diffusion=diffusion,
        )
        for space, aspect in geometries
        for decay in (0.42, 0.90)
        for damping in (0.25, 0.82)
        for diffusion in (0.18, 0.62)
    )


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
    controls: tuple[VelvetControls, ...],
    control_batch: int,
) -> tuple[float, float]:
    values = []
    worst_band = 0.0
    for first in range(0, len(controls), control_batch):
        wall, resolvent = model.response(
            frequencies, controls[first : first + control_batch]
        )
        combined = torch.cat((wall.flatten(2), resolvent.flatten(2)), dim=2)
        per_control, diagnostics = resonance_per_control(combined, frequencies)
        values.append(per_control)
        worst_band = max(worst_band, float(diagnostics["worst_band"]))
    return float(aggregate_resonance(torch.cat(values))), worst_band


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", nargs="?", type=Path)
    parser.add_argument("--seconds", type=float, default=8.0)
    parser.add_argument("--control-batch", type=int, default=5)
    parser.add_argument("--blends", type=float, nargs="+", default=(1.0,))
    parser.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu"
    )
    args = parser.parse_args()
    if args.seconds <= 0.0 or args.control_batch < 1:
        parser.error("seconds and control-batch must be positive")
    if any(blend <= 0.0 or blend > 1.0 for blend in args.blends):
        parser.error("blends must be greater than zero and at most one")

    device = torch.device(args.device)
    frequencies = torch.arange(20.0, 20_000.0, 1.0 / args.seconds, device=device)
    baseline = DifferentiableVelvetReverb().to(device=device)
    models = [("Base", baseline)]
    if args.artifact is not None:
        for blend in args.blends:
            optimized = DifferentiableVelvetReverb().to(device=device)
            load_artifact(optimized, args.artifact)
            with torch.no_grad():
                ratios = torch.lerp(
                    baseline.base_main_ratios,
                    optimized.base_main_ratios,
                    blend,
                )
                optimized.base_main_ratios.copy_(ratios / ratios.mean())
            name = "Optimized" if blend == 1.0 else f"Blend {blend:.3g}"
            models.append((name, optimized))
    grids = (
        ("training", control_grid()),
        ("held-out", validation_control_grid()),
        ("secondary-validation", secondary_validation_control_grid()),
    )
    for name, model in models:
        for grid_name, controls in grids:
            objective, worst_band = score(
                model, frequencies, controls, args.control_batch
            )
            print(
                f"{name} {grid_name}: controls={len(controls)} "
                f"objective={objective:.9g} worst-band={worst_band:.9g}"
            )


if __name__ == "__main__":
    main()
