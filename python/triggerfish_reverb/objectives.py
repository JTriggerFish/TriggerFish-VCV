"""Reference-free perceptual objectives for the velvet reverb optimizer."""

from __future__ import annotations

import math

import torch
import torch.nn.functional as F

from .velvet import smooth_worst


def _physical_spectral_residual(
    channels: torch.Tensor,
    frequencies: torch.Tensor,
    window_hz: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Remove a local spectral mean using a grid-invariant Hz window.

    Broad windows are evaluated on a proportionally decimated grid. This
    bounds work and storage while retaining roughly 64 samples per smoothing
    window at every scale.
    """
    step_hz = float((frequencies[1] - frequencies[0]).detach().cpu())
    stride = max(1, int(round(window_hz / (64.0 * step_hz))))
    sampled = channels[..., ::stride]
    sampled_frequencies = frequencies[::stride]
    kernel = max(3, int(round(window_hz / (step_hz * stride))))
    if kernel % 2 == 0:
        kernel += 1
    radius = kernel // 2
    local_mean = F.avg_pool1d(
        F.pad(sampled, (radius, radius), mode="replicate"),
        kernel_size=kernel,
        stride=1,
    )
    return sampled - local_mean, sampled_frequencies


def resonance_per_control(
    response: torch.Tensor, frequencies: torch.Tensor
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    """Penalize narrow modes, modal beating and broad comb coloration."""
    magnitude_db = (10.0 / math.log(10.0)) * torch.log(
        response.abs().square().clamp_min(1.0e-12)
    )
    channels = magnitude_db.flatten(2).transpose(1, 2)
    band_edges = (
        20.0,
        80.0,
        160.0,
        320.0,
        640.0,
        1_280.0,
        2_560.0,
        5_120.0,
        10_240.0,
        20_000.0,
    )
    band_losses = []
    # Metallic tails can arise from sub-Hz pole prominence, several-Hz modal
    # beating, or wider comb coloration. Measure every octave at all physical
    # scales so grid duration cannot change the definition of the objective.
    for window_hz, threshold_db, scale_weight in (
        (2.0, 1.2, 1.0),
        (8.0, 1.5, 1.0),
        (32.0, 2.0, 0.75),
        (128.0, 3.0, 0.40),
    ):
        residual, sampled_frequencies = _physical_spectral_residual(
            channels, frequencies, window_hz
        )
        for lower, upper in zip(band_edges[:-1], band_edges[1:], strict=True):
            selected = (sampled_frequencies >= lower) & (sampled_frequencies < upper)
            if not torch.any(selected):
                continue
            values = residual[..., selected]
            # Metallic ringing is caused by locally prominent poles. Deep
            # transfer zeros are often benign spatial cancellation and must
            # not consume the same loss budget as an upward resonance.
            excess = F.relu(values - threshold_db)
            rms = excess.square().mean(dim=(1, 2))
            peak = torch.logsumexp(0.35 * excess.flatten(1), dim=1) / 0.35
            peak = peak - math.log(excess[0].numel()) / 0.35
            band_losses.append(
                scale_weight * (rms + 0.12 * F.relu(peak - 2.5).square())
            )
    bands = torch.stack(band_losses, dim=1)
    per_control = bands.mean(dim=1) + 0.75 * torch.stack(
        [smooth_worst(item) for item in bands]
    )
    return per_control, {
        "resonance": per_control.mean(),
        "worst_band": bands.max(),
    }


def aggregate_resonance(per_control: torch.Tensor) -> torch.Tensor:
    """Combine controls only after every corner has been evaluated."""
    return per_control.mean() + 0.75 * smooth_worst(per_control)


def resonance_loss(
    response: torch.Tensor, frequencies: torch.Tensor
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    per_control, diagnostics = resonance_per_control(response, frequencies)
    return aggregate_resonance(per_control), diagnostics
