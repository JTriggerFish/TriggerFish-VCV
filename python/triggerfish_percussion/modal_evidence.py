"""Multiband, multiresolution corroboration of damped-mode candidates."""

from dataclasses import dataclass

import numpy as np
from scipy.signal import firwin, lfilter

from .modes import DampedMode, estimate_damped_modes


@dataclass(frozen=True)
class EspritPass:
    minimum_hz: float
    maximum_hz: float
    mode_count: int
    pencil_samples: int
    maximum_analysis_samples: int = 8192
    filter_taps: int = 257


@dataclass(frozen=True)
class ModeEvidence:
    mode: DampedMode
    observation_count: int
    hit_count: int
    frequency_std_cents: float
    decay_std_log: float
    estimates: tuple[DampedMode, ...]


@dataclass
class _Cluster:
    estimates: list[DampedMode]
    sources: set[tuple[int, int]]
    hits: set[int]


def estimate_modal_evidence(
    hits: tuple[np.ndarray, ...],
    sample_rate: float,
    passes: tuple[EspritPass, ...],
    merge_cents: float = 20.0,
) -> tuple[ModeEvidence, ...]:
    if not hits or not passes or merge_cents <= 0:
        raise ValueError("modal evidence requires hits, passes, and a merge tolerance")
    observations: list[tuple[DampedMode, tuple[int, int]]] = []
    for hit_index, hit in enumerate(hits):
        for pass_index, configuration in enumerate(passes):
            modes = estimate_subband_modes(hit, sample_rate, configuration)
            observations.extend((mode, (hit_index, pass_index)) for mode in modes)
    clusters = _cluster_observations(observations, merge_cents)
    evidence = tuple(_summarize(cluster) for cluster in clusters)
    return tuple(sorted(evidence, key=lambda item: item.mode.frequency_hz))


def estimate_subband_modes(
    samples: np.ndarray, sample_rate: float, configuration: EspritPass
) -> tuple[DampedMode, ...]:
    signal = np.asarray(samples, dtype=np.float64)
    if not 0 <= configuration.minimum_hz < configuration.maximum_hz < 0.5 * sample_rate:
        raise ValueError("ESPRIT pass band must lie below Nyquist")
    if configuration.filter_taps < 3 or configuration.filter_taps % 2 == 0:
        raise ValueError("ESPRIT subband filter needs an odd tap count")
    cutoff: float | list[float]
    pass_zero: bool
    if configuration.minimum_hz == 0:
        cutoff = configuration.maximum_hz
        pass_zero = True
    else:
        cutoff = [configuration.minimum_hz, configuration.maximum_hz]
        pass_zero = False
    coefficients = firwin(
        configuration.filter_taps, cutoff, pass_zero=pass_zero, fs=sample_rate
    )
    if signal.size <= configuration.filter_taps + 16:
        raise ValueError("ESPRIT hit is too short for the selected subband filter")
    filtered = lfilter(coefficients, [1.0], signal)
    filtered = filtered[configuration.filter_taps - 1 :]
    modes = estimate_damped_modes(
        filtered,
        sample_rate,
        configuration.mode_count,
        configuration.pencil_samples,
        configuration.maximum_analysis_samples,
    )
    return tuple(
        mode
        for mode in modes
        if configuration.minimum_hz <= mode.frequency_hz <= configuration.maximum_hz
    )


def _cluster_observations(
    observations: list[tuple[DampedMode, tuple[int, int]]], merge_cents: float
) -> list[_Cluster]:
    clusters: list[_Cluster] = []
    for mode, source in sorted(observations, key=lambda item: item[0].frequency_hz):
        available = [cluster for cluster in clusters if source not in cluster.sources]
        distances = [
            abs(_cents(mode.frequency_hz, _center(cluster))) for cluster in available
        ]
        if distances and min(distances) <= merge_cents:
            cluster = available[int(np.argmin(distances))]
            cluster.estimates.append(mode)
            cluster.sources.add(source)
            cluster.hits.add(source[0])
        else:
            clusters.append(_Cluster([mode], {source}, {source[0]}))
    return clusters


def _summarize(cluster: _Cluster) -> ModeEvidence:
    frequencies = np.array([mode.frequency_hz for mode in cluster.estimates])
    decays = np.array([mode.decay_seconds for mode in cluster.estimates])
    amplitudes = np.array([mode.amplitude for mode in cluster.estimates])
    frequency = float(np.exp(np.mean(np.log(frequencies))))
    decay = float(np.exp(np.mean(np.log(decays))))
    phase = float(cluster.estimates[0].phase_radians)
    frequency_spread = float(np.std(1200.0 * np.log2(frequencies / frequency)))
    decay_spread = float(np.std(np.log(decays / decay)))
    representative = DampedMode(frequency, decay, float(np.median(amplitudes)), phase)
    return ModeEvidence(
        representative,
        len(cluster.estimates),
        len(cluster.hits),
        frequency_spread,
        decay_spread,
        tuple(cluster.estimates),
    )


def _center(cluster: _Cluster) -> float:
    return float(
        np.exp(np.mean(np.log([mode.frequency_hz for mode in cluster.estimates])))
    )


def _cents(first_hz: float, second_hz: float) -> float:
    return float(1200.0 * np.log2(first_hz / second_hz))
