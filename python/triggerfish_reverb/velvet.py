"""Differentiable frequency-domain model of the current production VFM tail.

The model mirrors the static, unmodulated and unshifted C++ late-field loop.
Runtime crossfades, modulation and shimmer are intentionally excluded because
they do not define the time-invariant room poles optimized here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import torch
from torch import Tensor, nn

LINE_COUNT = 16
WALL_COUNT = 6
VELVET_STAGE_COUNT = 2
SPEED_OF_SOUND = 343.0
REFERENCE_ROOM_DIMENSIONS_METRES = (7.09929574, 9.35414347, 4.38178046)

MAIN_DELAY_RATIO = (
    0.5668773,
    0.6296269,
    0.6799816,
    0.7295616,
    0.7899872,
    0.8504128,
    0.9201346,
    0.9681652,
    1.0258794,
    1.0773961,
    1.1362723,
    1.2001840,
    1.2629336,
    1.3280074,
    1.3864962,
    1.4480839,
)
VELVET_DELAY_MS = (
    (
        0.125,
        0.875,
        1.375,
        2.0,
        2.75,
        3.25,
        3.875,
        4.5,
        5.25,
        5.75,
        6.375,
        7.125,
        7.75,
        8.25,
        9.125,
        9.75,
    ),
    (
        0.375,
        1.0,
        1.625,
        2.375,
        2.875,
        3.625,
        4.125,
        4.875,
        5.375,
        6.125,
        6.75,
        7.375,
        8.0,
        8.75,
        9.375,
        10.125,
    ),
)
PERMUTATIONS = (
    (11, 1, 15, 0, 6, 3, 4, 13, 8, 10, 14, 12, 9, 7, 2, 5),
    (9, 10, 12, 4, 7, 5, 13, 1, 3, 6, 0, 8, 2, 15, 11, 14),
    (5, 12, 1, 9, 14, 3, 10, 0, 7, 15, 4, 11, 2, 13, 8, 6),
)
SIGNS = (
    (-1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1),
    (1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1),
    (1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1),
)


@dataclass(frozen=True)
class VelvetControls:
    """Normalized controls that affect the static production late-field loop."""

    space: float = 0.5
    aspect: float = 0.5
    decay: float = 0.55
    damping: float = 0.18
    diffusion: float = 0.75


def _walsh_matrix(lines: int = LINE_COUNT) -> Tensor:
    result = torch.empty(lines, lines)
    for row in range(lines):
        for column in range(lines):
            result[row, column] = -1.0 if (row & column).bit_count() & 1 else 1.0
    return result / math.sqrt(lines)


def _wall_projection(lines: int = LINE_COUNT, walls: int = WALL_COUNT) -> Tensor:
    return _walsh_matrix(lines)[:, :walls]


def _smoothstep(value: Tensor) -> Tensor:
    limited = value.clamp(0.0, 1.0)
    return limited.square() * (3.0 - 2.0 * limited)


def room_dimensions_metres(space: Tensor, aspect: Tensor) -> Tensor:
    """Mirror MakeEarlyReflectionRoom() and RoomReverb::MakeRoom()."""
    minimum = space.new_tensor((2.8, 3.5, 2.4))
    maximum = space.new_tensor((18.0, 25.0, 8.0))
    amount = _smoothstep(space).unsqueeze(-1)
    dimensions = torch.exp(
        torch.log(minimum) + amount * (torch.log(maximum) - torch.log(minimum))
    )
    aspect_ratio = torch.exp((2.0 * aspect.clamp(0.0, 1.0) - 1.0) * math.log(1.8))
    root_aspect = torch.sqrt(aspect_ratio)
    dimensions = dimensions.clone()
    dimensions[:, 0] *= root_aspect
    dimensions[:, 1] /= root_aspect
    return dimensions


def mean_free_time_seconds(dimensions: Tensor) -> Tensor:
    """Mirror LateReverb::MeanFreeTimeSeconds()."""
    limited = dimensions.clamp(1.0, 40.0)
    length, width, height = limited.unbind(dim=-1)
    volume = length * width * height
    surface = 2.0 * (length * width + length * height + width * height)
    return (4.0 * volume / (surface * SPEED_OF_SOUND)).clamp(0.003, 0.078)


class DifferentiableVelvetReverb(nn.Module):
    """Current two-stage production loop with trainable delay ratios."""

    def __init__(self, sample_rate: float = 48_000.0):
        super().__init__()
        if not math.isfinite(sample_rate) or sample_rate <= 0.0:
            raise ValueError("sample_rate must be positive and finite")
        self.sample_rate = float(sample_rate)
        self.register_buffer("base_main_ratios", torch.tensor(MAIN_DELAY_RATIO))
        self.register_buffer("base_velvet_ms", torch.tensor(VELVET_DELAY_MS))
        self.register_buffer("permutations", torch.tensor(PERMUTATIONS))
        self.register_buffer("signs", torch.tensor(SIGNS, dtype=torch.float32))
        self.register_buffer("hadamard", _walsh_matrix())
        self.register_buffer("wall_projection", _wall_projection())
        # Fine-tune around the proven prime-delay heuristic without allowing
        # adjacent main paths to collapse into a slow beating pair.
        self.register_buffer("main_ratio_search_radius", torch.tensor(0.005))
        self.raw_main_ratios = nn.Parameter(torch.zeros(LINE_COUNT))

    @property
    def complex_dtype(self) -> torch.dtype:
        return (
            torch.complex128
            if self.raw_main_ratios.dtype == torch.float64
            else torch.complex64
        )

    def main_delay_ratios(self) -> Tensor:
        ratios = self.base_main_ratios + self.main_ratio_search_radius * torch.tanh(
            self.raw_main_ratios
        )
        ratios = torch.sort(ratios).values
        return ratios / ratios.mean()

    def coefficient_tensors(self) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        return (
            self.main_delay_ratios(),
            self.base_velvet_ms,
            self.permutations,
            self.signs,
        )

    @torch.no_grad()
    def set_coefficients(
        self,
        main_delay_ratio: Tensor,
        velvet_delay_ms: Tensor,
        permutations: Tensor,
        signs: Tensor,
    ) -> None:
        if main_delay_ratio.shape != (LINE_COUNT,):
            raise ValueError("main_delay_ratio must contain 16 values")
        if velvet_delay_ms.shape != (VELVET_STAGE_COUNT, LINE_COUNT):
            raise ValueError("velvet_delay_ms must have shape (2, 16)")
        if permutations.shape != (VELVET_STAGE_COUNT + 1, LINE_COUNT):
            raise ValueError("permutations must have shape (3, 16)")
        if signs.shape != (VELVET_STAGE_COUNT + 1, LINE_COUNT):
            raise ValueError("signs must have shape (3, 16)")
        self.base_main_ratios.copy_(main_delay_ratio / main_delay_ratio.mean())
        self.base_velvet_ms.copy_(velvet_delay_ms)
        self.permutations.copy_(permutations)
        self.signs.copy_(signs)
        self.raw_main_ratios.zero_()

    def transform(self, index: int) -> Tensor:
        """Return the exact fixed signed-Hadamard transform used by C++."""
        if not 0 <= index <= VELVET_STAGE_COUNT:
            raise IndexError("transform index is outside the current topology")
        return self.signs[index].unsqueeze(-1) * self.hadamard[self.permutations[index]]

    def _control_tensors(
        self, controls: tuple[VelvetControls, ...]
    ) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        if not controls:
            raise ValueError("at least one control point is required")
        reference = self.raw_main_ratios
        return tuple(
            reference.new_tensor([getattr(item, name) for item in controls])
            for name in ("space", "aspect", "decay", "damping", "diffusion")
        )

    def room_state(
        self, controls: tuple[VelvetControls, ...]
    ) -> tuple[Tensor, Tensor, Tensor]:
        space, aspect, _, _, diffusion = self._control_tensors(controls)
        dimensions = room_dimensions_metres(space, aspect)
        mean_time = mean_free_time_seconds(dimensions)
        reference_dimensions = space.new_tensor(
            REFERENCE_ROOM_DIMENSIONS_METRES
        ).reshape(1, 3)
        reference_time = mean_free_time_seconds(reference_dimensions)[0]
        room_scale = (mean_time / reference_time).clamp(0.35, 2.25)
        diffusion_scale = 0.20 + 1.30 * _smoothstep(diffusion)
        return mean_time, room_scale, diffusion_scale

    def velvet_delay_samples(self, controls: tuple[VelvetControls, ...]) -> Tensor:
        _, room_scale, diffusion_scale = self.room_state(controls)
        samples = (
            self.base_velvet_ms.reshape(1, VELVET_STAGE_COUNT, LINE_COUNT)
            * room_scale.reshape(-1, 1, 1)
            * diffusion_scale.reshape(-1, 1, 1)
            * self.sample_rate
            / 1_000.0
        )
        return torch.round(samples).clamp_min(1.0)

    def _one_pole_lowpass(self, frequencies_hz: Tensor, cutoff_hz: float) -> Tensor:
        alpha = 1.0 - math.exp(-2.0 * math.pi * cutoff_hz / self.sample_rate)
        omega = 2.0 * math.pi * frequencies_hz / self.sample_rate
        delay = torch.exp(-1j * omega).to(self.complex_dtype)
        return alpha / (1.0 - (1.0 - alpha) * delay)

    def _decay_times(
        self, controls: tuple[VelvetControls, ...]
    ) -> tuple[Tensor, Tensor, Tensor]:
        _, _, decay, damping, _ = self._control_tensors(controls)
        damping = _smoothstep(damping)
        mid = 0.25 * torch.exp(decay.clamp(0.0, 1.0) * math.log(32.0))
        high = mid * torch.pow(mid.new_tensor(0.20), damping)
        low = mid * torch.pow(mid.new_tensor(0.72), damping)
        return low, mid, high

    def multiband_loss_response(
        self,
        frequencies_hz: Tensor,
        path_seconds: Tensor,
        controls: tuple[VelvetControls, ...],
    ) -> Tensor:
        """Exact transfer of MultibandDecayFilter for each path."""
        low_t60, mid_t60, high_t60 = self._decay_times(controls)

        def gain(t60: Tensor) -> Tensor:
            return torch.pow(
                path_seconds.new_tensor(10.0),
                -3.0 * path_seconds / t60.reshape(-1, 1).clamp_min(1.0e-4),
            )

        low_gain = gain(low_t60).unsqueeze(1)
        mid_gain = gain(mid_t60).unsqueeze(1)
        high_gain = gain(high_t60).unsqueeze(1)
        lowpass = self._one_pole_lowpass(frequencies_hz, 220.0).reshape(1, -1, 1)
        below_high = self._one_pole_lowpass(
            frequencies_hz, min(3_500.0, 0.20 * self.sample_rate)
        ).reshape(1, -1, 1)
        return (
            low_gain * lowpass
            + mid_gain * (below_high - lowpass)
            + high_gain * (1.0 - below_high)
        )

    def _fractional_delay_response(
        self, frequencies_hz: Tensor, delay_samples: Tensor
    ) -> Tensor:
        """Frequency response of the runtime four-point Lagrange read head."""
        integer = torch.floor(delay_samples)
        mu = delay_samples - integer
        coefficients = torch.stack(
            (
                -mu * (mu - 1.0) * (mu - 2.0) / 6.0,
                (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0,
                -(mu + 1.0) * mu * (mu - 2.0) / 2.0,
                (mu + 1.0) * mu * (mu - 1.0) / 6.0,
            ),
            dim=-1,
        )
        offsets = delay_samples.new_tensor((-1.0, 0.0, 1.0, 2.0))
        distances = integer.unsqueeze(-1) + offsets
        omega = 2.0 * math.pi * frequencies_hz / self.sample_rate
        phase = torch.exp(-1j * omega.reshape(1, -1, 1, 1) * distances.unsqueeze(1)).to(
            self.complex_dtype
        )
        return (coefficients.unsqueeze(1) * phase).sum(dim=-1)

    def feedback_response(
        self,
        frequencies_hz: Tensor,
        controls: tuple[VelvetControls, ...],
        *,
        include_losses: bool = True,
    ) -> Tensor:
        """Return U2 L2 D2 U1 L1 D1 U0 for every control and frequency."""
        frequencies_hz = frequencies_hz.to(self.raw_main_ratios)
        batch = len(controls)
        frequency_count = frequencies_hz.numel()
        result = (
            self.transform(0)
            .to(self.complex_dtype)
            .reshape(1, 1, LINE_COUNT, LINE_COUNT)
        )
        result = result.expand(batch, frequency_count, -1, -1)
        stage_samples = self.velvet_delay_samples(controls)
        omega = 2.0 * math.pi * frequencies_hz / self.sample_rate
        for stage in range(VELVET_STAGE_COUNT):
            samples = stage_samples[:, stage, :]
            diagonal = torch.exp(
                -1j * omega.reshape(1, -1, 1) * samples.unsqueeze(1)
            ).to(self.complex_dtype)
            if include_losses:
                diagonal = diagonal * self.multiband_loss_response(
                    frequencies_hz, samples / self.sample_rate, controls
                )
            result = diagonal.unsqueeze(-1) * result
            transform = self.transform(stage + 1).to(self.complex_dtype)
            result = torch.einsum("ij,bfjk->bfik", transform, result)
        return result

    def response(
        self,
        frequencies_hz: Tensor,
        controls: tuple[VelvetControls, ...],
    ) -> tuple[Tensor, Tensor]:
        """Return six-wall transfer and complete internal loop resolvent."""
        if frequencies_hz.ndim != 1:
            raise ValueError("frequencies_hz must be one-dimensional")
        frequencies_hz = frequencies_hz.to(self.raw_main_ratios)
        mean_time, _, _ = self.room_state(controls)
        ratios = self.main_delay_ratios()
        main_samples = (
            mean_time.reshape(-1, 1) * self.sample_rate * ratios.reshape(1, -1)
        )
        main_delay = self._fractional_delay_response(frequencies_hz, main_samples)
        main_loss = self.multiband_loss_response(
            frequencies_hz, main_samples / self.sample_rate, controls
        )
        feedback = self.feedback_response(frequencies_hz, controls, include_losses=True)
        feedback = feedback * main_loss.unsqueeze(-2)
        loop = main_delay.unsqueeze(-1) * feedback
        identity = torch.eye(LINE_COUNT, dtype=self.complex_dtype, device=loop.device)
        system = identity - loop
        resolvent = torch.linalg.solve(system, identity)
        injection = 0.42 * self.wall_projection.to(self.complex_dtype)
        right = main_delay.unsqueeze(-1) * injection.reshape(
            1, 1, LINE_COUNT, WALL_COUNT
        )
        internal = torch.einsum("bfij,bfjs->bfis", resolvent, right)
        wall_output = torch.einsum(
            "lw,bfls->bfws", self.wall_projection.to(internal.dtype), internal
        )
        return wall_output, resolvent


def smooth_worst(values: Tensor, sharpness: float = 12.0) -> Tensor:
    return (
        torch.logsumexp(sharpness * values, dim=0) - math.log(values.numel())
    ) / sharpness
