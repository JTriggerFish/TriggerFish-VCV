"""Report-neutral reference/synthesis comparison contract."""

from dataclasses import dataclass

import numpy as np

from .alignment import Alignment, measure_alignment, shift_with_zeros
from .audio_io import AudioBuffer
from .distances import (
    LossTerm,
    energy_matching_gain,
    erb_trajectory_distance,
    log_spectral_distance,
)
from .erb import ErbFilterbank
from .transforms import StftConfig, StftResult, stft


@dataclass(frozen=True)
class ComparisonConfig:
    stft_configs: tuple[StftConfig, ...]
    erb_minimum_hz: float = 40.0
    erb_maximum_hz: float = 18_000.0
    erb_band_count: int = 40
    level_mode: str = "preserve"
    minimum_region_overlap: float = 1.0
    regions: tuple[tuple[str, float, float | None], ...] = (
        ("contact", 0.0, 0.015),
        ("bloom", 0.015, 0.120),
        ("early_body", 0.120, 0.600),
        ("tail", 0.600, None),
    )


@dataclass(frozen=True)
class SpectralComparison:
    reference: StftResult
    candidate: StftResult
    erb_centers_hz: np.ndarray
    reference_erb_power: np.ndarray
    candidate_erb_power: np.ndarray
    losses: tuple[LossTerm, ...]


@dataclass(frozen=True)
class Comparison:
    reference: AudioBuffer
    candidate: AudioBuffer
    alignment: Alignment
    candidate_gain: float
    resolutions: tuple[SpectralComparison, ...]

    @property
    def losses(self) -> tuple[LossTerm, ...]:
        return tuple(
            loss for resolution in self.resolutions for loss in resolution.losses
        )


def compare_audio(
    reference: AudioBuffer,
    candidate: AudioBuffer,
    config: ComparisonConfig,
) -> Comparison:
    ref = reference.mono()
    current = candidate.mono()
    if ref.sample_rate != current.sample_rate:
        raise ValueError("comparison audio must share one sample rate")
    if not config.stft_configs:
        raise ValueError("comparison requires at least one STFT resolution")
    if not 0.0 < config.minimum_region_overlap <= 1.0:
        raise ValueError("minimum region overlap must lie in (0, 1]")
    alignment = measure_alignment(ref.samples, current.samples, ref.sample_rate)
    shifted = shift_with_zeros(
        _same_length(current.samples, ref.samples.size), alignment.candidate_lag_samples
    )
    gain = _comparison_gain(ref.samples, shifted, config.level_mode)
    aligned = AudioBuffer(gain * shifted, ref.sample_rate)
    onset_seconds = alignment.reference_onset / ref.sample_rate
    resolutions = tuple(
        _compare_resolution(ref, aligned, stft_config, config, onset_seconds)
        for stft_config in config.stft_configs
    )
    return Comparison(ref, aligned, alignment, gain, resolutions)


def _compare_resolution(
    reference: AudioBuffer,
    candidate: AudioBuffer,
    stft_config: StftConfig,
    config: ComparisonConfig,
    onset_seconds: float,
) -> SpectralComparison:
    ref_stft = stft(reference.samples, reference.sample_rate, stft_config)
    candidate_stft = stft(candidate.samples, candidate.sample_rate, stft_config)
    maximum_hz = min(config.erb_maximum_hz, float(ref_stft.frequencies_hz[-1]))
    bank = ErbFilterbank.create(
        ref_stft.frequencies_hz,
        config.erb_minimum_hz,
        maximum_hz,
        config.erb_band_count,
    )
    reference_erb = bank.apply_power(ref_stft.power)
    candidate_erb = bank.apply_power(candidate_stft.power)
    losses: list[LossTerm] = []
    resolution = f"w{stft_config.window_samples}"
    for region, start, end in config.regions:
        selected = _region_frames(
            ref_stft.times_seconds,
            stft_config.window_samples / reference.sample_rate,
            onset_seconds + start,
            None if end is None else onset_seconds + end,
            reference.samples.size / reference.sample_rate,
            config.minimum_region_overlap,
        )
        if not np.any(selected):
            continue
        spectral = log_spectral_distance(
            ref_stft.magnitude[:, selected], candidate_stft.magnitude[:, selected]
        )
        erb_losses = erb_trajectory_distance(
            reference_erb[:, selected], candidate_erb[:, selected]
        )
        losses.extend(
            _scoped(loss, region, resolution) for loss in (spectral, *erb_losses)
        )
    return SpectralComparison(
        ref_stft,
        candidate_stft,
        bank.centers_hz,
        reference_erb,
        candidate_erb,
        tuple(losses),
    )


def _scoped(loss: LossTerm, region: str, resolution: str) -> LossTerm:
    return LossTerm(f"{loss.name}/{region}/{resolution}", loss.value, loss.unit)


def _comparison_gain(reference: np.ndarray, candidate: np.ndarray, mode: str) -> float:
    if mode == "preserve":
        return 1.0
    if mode == "energy":
        return energy_matching_gain(reference, candidate)
    raise ValueError(f"unknown comparison level mode: {mode}")


def _same_length(samples: np.ndarray, length: int) -> np.ndarray:
    result = np.zeros(length, dtype=np.float64)
    count = min(length, samples.size)
    result[:count] = samples[:count]
    return result


def _region_frames(
    frame_centers: np.ndarray,
    window_seconds: float,
    region_start: float,
    region_end: float | None,
    signal_end: float,
    minimum_overlap: float,
) -> np.ndarray:
    half_window = 0.5 * window_seconds
    frame_start = frame_centers - half_window
    frame_end = frame_centers + half_window
    bounded_end = signal_end if region_end is None else min(region_end, signal_end)
    overlap = np.maximum(
        0.0,
        np.minimum(frame_end, bounded_end) - np.maximum(frame_start, region_start),
    )
    fraction = overlap / window_seconds
    return fraction >= minimum_overlap - 1.0e-12
