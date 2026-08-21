#!/usr/bin/env python3
"""Optimize coefficients for the current two-stage production VFM.

The objective is reference-free. Every accepted continuous step is evaluated
over the full audible grid and the complete declared control grid. Control
sub-batches change peak memory only; they preserve the exact aggregate gradient.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import torch
import torch.nn.functional as F

from triggerfish_reverb.objectives import (
    aggregate_resonance,
    resonance_per_control,
)
from triggerfish_reverb.velvet import (
    LINE_COUNT,
    VELVET_STAGE_COUNT,
    DifferentiableVelvetReverb,
    VelvetControls,
)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_FORMAT = "TriggerFish current two-stage VFM coefficients v1"
DEFAULT_OUTPUT = ROOT / "data/reverb-calibration/velvet-vfm-current-v1.json"


def control_grid() -> tuple[VelvetControls, ...]:
    """Cover all loss/diffusion corners at five distinct room geometries."""
    geometries = (
        (0.0, 0.5),
        (0.5, 0.0),
        (0.5, 0.5),
        (0.5, 1.0),
        (1.0, 0.5),
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
        for decay in (0.0, 0.55, 1.0)
        for damping in (0.0, 0.18, 1.0)
        for diffusion in (0.0, 0.75, 1.0)
    )


def model_resonance_per_control(
    model: DifferentiableVelvetReverb,
    frequencies: torch.Tensor,
    controls: tuple[VelvetControls, ...],
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    wall_response, resolvent = model.response(frequencies, controls)
    # The wall transfer measures the physically reachable paths. The complete
    # resolvent prevents a pole outside the default six-wall subspace from
    # escaping the objective.
    combined = torch.cat((wall_response.flatten(2), resolvent.flatten(2)), dim=2)
    return resonance_per_control(combined, frequencies)


@torch.no_grad()
def evaluate_resonance_batches(
    model: DifferentiableVelvetReverb,
    frequencies: torch.Tensor,
    controls: tuple[VelvetControls, ...],
    batch_size: int,
) -> tuple[torch.Tensor, float]:
    values = []
    worst_band = 0.0
    for first in range(0, len(controls), batch_size):
        batch = controls[first : first + batch_size]
        per_control, diagnostics = model_resonance_per_control(
            model, frequencies, batch
        )
        values.append(per_control)
        worst_band = max(worst_band, float(diagnostics["worst_band"]))
    return torch.cat(values), worst_band


def global_control_weights(per_control: torch.Tensor) -> torch.Tensor:
    """Exact derivative weights of mean + 0.75 * smooth-worst."""
    return torch.full_like(
        per_control, 1.0 / per_control.numel()
    ) + 0.75 * torch.softmax(12.0 * per_control, dim=0)


def delay_regularization(model: DifferentiableVelvetReverb) -> torch.Tensor:
    """Discourage coincident delays and repeated path-length differences."""
    ratios = model.main_delay_ratios()
    gaps = torch.diff(ratios)
    separation = F.relu(0.035 - gaps).square().mean() / (0.035**2)
    pairwise = (ratios[:, None] - ratios[None, :]).abs()
    pairwise = pairwise[
        torch.triu(torch.ones_like(pairwise, dtype=torch.bool), diagonal=1)
    ]
    repeated = (pairwise[:, None] - pairwise[None, :]).abs()
    repeated = repeated[
        torch.triu(torch.ones_like(repeated, dtype=torch.bool), diagonal=1)
    ]
    repeated_loss = torch.exp(-(repeated / 0.0035).square()).mean()

    velvet_samples = model.base_velvet_ms * model.sample_rate / 1_000.0
    velvet_separation = sum(
        F.relu(3.0 - torch.diff(torch.sort(stage).values)).square().mean() / 9.0
        for stage in velvet_samples
    )
    return separation + 5.0 * repeated_loss + 0.1 * velvet_separation


def complete_objective(
    model: DifferentiableVelvetReverb,
    per_control: torch.Tensor,
) -> torch.Tensor:
    return aggregate_resonance(per_control) + 0.12 * delay_regularization(model)


@torch.no_grad()
def screen_discrete_coefficients(
    model: DifferentiableVelvetReverb,
    candidates: int,
    seed: int,
    control_batch: int,
) -> tuple[int, float]:
    """Screen signed permutations and integer-delay layouts."""
    if candidates <= 1:
        return seed, float("nan")
    device = model.raw_main_ratios.device
    generator = torch.Generator(device=device).manual_seed(seed)
    frequencies = torch.arange(
        20.0,
        20_000.0,
        1.0,
        device=device,
        dtype=model.raw_main_ratios.dtype,
    )
    controls = control_grid()
    initial_main = model.base_main_ratios.clone()
    initial_velvet = model.base_velvet_ms.clone()
    best_main = initial_main.clone()
    best_velvet = initial_velvet.clone()
    best_permutations = model.permutations.clone()
    best_signs = model.signs.clone()
    best_seed = seed
    best_loss = float("inf")

    for candidate in range(candidates):
        candidate_seed = seed + candidate
        if candidate > 0:
            for transform in range(VELVET_STAGE_COUNT + 1):
                model.permutations[transform] = torch.randperm(
                    LINE_COUNT, generator=generator, device=device
                )
                bits = torch.randint(
                    0,
                    2,
                    (LINE_COUNT,),
                    generator=generator,
                    device=device,
                    dtype=torch.int64,
                )
                model.signs[transform] = 1.0 - 2.0 * bits.to(model.signs.dtype)

            jitter = 0.025 * (
                2.0
                * torch.rand(
                    LINE_COUNT,
                    generator=generator,
                    device=device,
                    dtype=initial_main.dtype,
                )
                - 1.0
            )
            candidate_main = torch.sort(initial_main + jitter).values
            model.base_main_ratios.copy_(candidate_main / candidate_main.mean())

            base_samples = torch.round(initial_velvet * model.sample_rate / 1_000.0)
            velvet_jitter = torch.randint(
                -6,
                7,
                base_samples.shape,
                generator=generator,
                device=device,
            )
            candidate_samples = torch.sort(
                (base_samples + velvet_jitter).clamp_min(1.0), dim=1
            ).values
            model.base_velvet_ms.copy_(candidate_samples * 1_000.0 / model.sample_rate)

        per_control, _ = evaluate_resonance_batches(
            model, frequencies, controls, control_batch
        )
        value = float(complete_objective(model, per_control))
        if value < best_loss:
            best_seed = candidate_seed
            best_loss = value
            best_main.copy_(model.base_main_ratios)
            best_velvet.copy_(model.base_velvet_ms)
            best_permutations.copy_(model.permutations)
            best_signs.copy_(model.signs)
        if candidate == 0 or (candidate + 1) % 8 == 0:
            print(
                f"screen={candidate + 1}/{candidates} "
                f"candidate={value:.6g} best={best_loss:.6g}",
                flush=True,
            )

    model.base_main_ratios.copy_(best_main)
    model.base_velvet_ms.copy_(best_velvet)
    model.permutations.copy_(best_permutations)
    model.signs.copy_(best_signs)
    model.raw_main_ratios.zero_()
    return best_seed, best_loss


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--seconds", type=float, default=32.0)
    parser.add_argument("--lr", type=float, default=0.20)
    parser.add_argument("--control-batch", type=int, default=5)
    parser.add_argument("--candidates", type=int, default=128)
    parser.add_argument("--seed", type=int, default=73021)
    parser.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu"
    )
    parser.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--diagnose-gradient", action="store_true")
    args = parser.parse_args()
    if (
        args.steps < 1
        or args.seconds <= 0.0
        or args.control_batch < 1
        or args.candidates < 1
        or args.lr <= 0.0
    ):
        parser.error(
            "steps, seconds, lr, control-batch and candidates must be positive"
        )

    torch.manual_seed(args.seed)
    device = torch.device(args.device)
    dtype = torch.float64 if args.dtype == "float64" else torch.float32
    model = DifferentiableVelvetReverb().to(device=device, dtype=dtype)
    selected_seed, screened_loss = screen_discrete_coefficients(
        model, args.candidates, args.seed, args.control_batch
    )
    print(
        f"selected discrete seed={selected_seed} screen_loss={screened_loss:.6g}",
        flush=True,
    )

    frequencies = torch.arange(
        20.0,
        20_000.0,
        1.0 / args.seconds,
        device=device,
        dtype=dtype,
    )
    controls = control_grid()

    @torch.no_grad()
    def evaluate_value() -> tuple[float, float]:
        per_control, worst_band = evaluate_resonance_batches(
            model, frequencies, controls, args.control_batch
        )
        return float(complete_objective(model, per_control)), worst_band

    def calculate_gradient() -> tuple[float, float, float]:
        model.raw_main_ratios.grad = None
        reference, worst_band = evaluate_resonance_batches(
            model, frequencies, controls, args.control_batch
        )
        weights = global_control_weights(reference).detach()
        for first in range(0, len(controls), args.control_batch):
            batch = controls[first : first + args.control_batch]
            per_control, _ = model_resonance_per_control(model, frequencies, batch)
            last = first + len(batch)
            (weights[first:last] * per_control).sum().backward()
        regularization = delay_regularization(model)
        (0.12 * regularization).backward()
        total = float(aggregate_resonance(reference))
        total += 0.12 * float(regularization.detach())
        return total, worst_band, float(regularization.detach())

    if args.diagnose_gradient:
        baseline, _, _ = calculate_gradient()
        gradient = model.raw_main_ratios.grad
        assert gradient is not None
        norm = gradient.norm().clamp_min(1.0e-20)
        original = model.raw_main_ratios.detach().clone()
        print(f"gradient baseline={baseline:.9g} norm={float(norm):.9g}")
        for distance in (1.0e-4, 1.0e-3, 1.0e-2, 0.05, 0.10):
            with torch.no_grad():
                model.raw_main_ratios.copy_(original - distance * gradient / norm)
            value, _ = evaluate_value()
            print(f"gradient distance={distance:.4g} objective={value:.9g}")
        with torch.no_grad():
            model.raw_main_ratios.copy_(original)
        return

    history: list[dict[str, float]] = []
    best_loss = float("inf")
    best_step = -1
    best_parameters: torch.Tensor | None = None

    for step in range(args.steps):
        total_value, worst_band, regularization = calculate_gradient()
        if total_value < best_loss:
            best_loss = total_value
            best_step = step
            best_parameters = model.raw_main_ratios.detach().cpu().clone()
        gradient = model.raw_main_ratios.grad
        assert gradient is not None and torch.isfinite(gradient).all()
        direction = gradient / gradient.norm().clamp_min(1.0e-20)
        original = model.raw_main_ratios.detach().clone()
        scheduled_radius = args.lr * (
            0.1 + 0.9 * 0.5 * (1.0 + math.cos(math.pi * step / args.steps))
        )
        accepted_distance = 0.0
        candidate_value = total_value
        for trial in range(8):
            distance = scheduled_radius * (0.5**trial)
            with torch.no_grad():
                model.raw_main_ratios.copy_(original - distance * direction)
            candidate_value, _ = evaluate_value()
            if candidate_value < total_value:
                accepted_distance = distance
                break
        if accepted_distance == 0.0:
            with torch.no_grad():
                model.raw_main_ratios.copy_(original)
            candidate_value = total_value

        history.append(
            {
                "step": step,
                "total": total_value,
                "accepted_total": candidate_value,
                "accepted_distance": accepted_distance,
                "worst_band": worst_band,
                "regularization": regularization,
            }
        )
        if step == 0 or (step + 1) % 10 == 0:
            memory = (
                torch.cuda.max_memory_allocated(device) / 2**30
                if device.type == "cuda"
                else 0.0
            )
            print(
                f"step={step + 1}/{args.steps} total={total_value:.6g} "
                f"accepted={candidate_value:.6g} trust={accepted_distance:.4g} "
                f"worst={worst_band:.6g} vram_gib={memory:.2f}",
                flush=True,
            )

    final_value, _ = evaluate_value()
    if final_value < best_loss:
        best_loss = final_value
        best_step = args.steps
        best_parameters = model.raw_main_ratios.detach().cpu().clone()
    assert best_parameters is not None
    with torch.no_grad():
        model.raw_main_ratios.copy_(best_parameters.to(device))

    main_ratios, velvet_ms, permutations, signs = model.coefficient_tensors()
    artifact = {
        "format": ARTIFACT_FORMAT,
        "architecture": {
            "line_count": LINE_COUNT,
            "velvet_stage_count": VELVET_STAGE_COUNT,
            "main_delay_scale": "room mean free time",
            "diffusion_law": "0.20 + 1.30 * smoothstep(control)",
            "loss": "three-band attenuation on every delay segment",
        },
        "sample_rate": model.sample_rate,
        "frequency_seconds": args.seconds,
        "frequency_range_hz": [20.0, 20_000.0],
        "control_grid_size": len(controls),
        "optimizer": {
            "name": "normalized steepest descent with backtracking",
            "initial_trust_radius": args.lr,
            "steps": args.steps,
        },
        "best_step": best_step,
        "best_loss": best_loss,
        "selected_discrete_seed": selected_seed,
        "main_delay_ratio": main_ratios.detach().cpu().tolist(),
        "velvet_delay_ms": velvet_ms.detach().cpu().tolist(),
        "permutations": permutations.detach().cpu().tolist(),
        "signs": signs.detach().cpu().tolist(),
        "history": history,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
