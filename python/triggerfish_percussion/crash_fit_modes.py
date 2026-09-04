"""Persistent-mode extraction, resynthesis, and sparse-bank fitting."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.optimize import nnls

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer
from .crash_fit_common import (
    CrashFitCell,
    fixed_rate,
    modal_power_features,
    render_sparse,
)
from .crash_model import CrashFit
from .modal_evidence import EspritPass, estimate_modal_evidence
from .modes import DampedMode, refit_mode_amplitudes, resynthesize_modes
from .t60_envelope import interpolate_t60


def fit_sparse_modes(
    cells: tuple[CrashFitCell, ...], initial: CrashFit
) -> tuple[CrashFit, dict[str, object]]:
    passes = modal_passes()
    tails = tuple(modal_tail(cell.reference) for cell in cells)
    required_hits = max(1, int(np.ceil(0.75 * len(cells))))
    stable = stable_modes(tails, 48000, passes, required_hits)
    ranked = stratified_modes(stable, len(initial.sparse_frequency_hz))
    selected = tuple(sorted(ranked, key=lambda mode: mode.frequency_hz))
    count = len(initial.sparse_frequency_hz)
    frequencies = list(initial.sparse_frequency_hz)
    amplitudes = [0.0] * count
    for index, mode in enumerate(selected):
        frequencies[index] = mode.frequency_hz
        amplitudes[index] = mode.amplitude
    peak = max(amplitudes, default=0.0)
    if peak > 0:
        amplitudes = [amplitude / peak for amplitude in amplitudes]
    fitted = replace(
        initial,
        sparse_frequency_hz=tuple(float(value) for value in frequencies),
        sparse_amplitude=tuple(float(value) for value in amplitudes),
        sparse_tune=1.0,
    )
    return fitted, {
        "required_hit_count": required_hits,
        "available_stable_mode_count": len(stable),
        "selected_mode_count": len(selected),
        "frequencies_hz": [mode.frequency_hz for mode in selected],
        "efold_decay_seconds": [mode.decay_seconds for mode in selected],
        "modal_t60_seconds": [np.log(1000.0) * mode.decay_seconds for mode in selected],
    }


def fit_sparse_projection(
    cells: tuple[CrashFitCell, ...], initial: CrashFit
) -> CrashFit:
    # This is only a linear anchor-energy initializer. Disable packet
    # turbulence and neighbour exchange while constructing basis columns:
    # otherwise zeroing all but one anchor changes the active state topology,
    # and those columns cannot be added to predict the full renderer.
    projection_model = replace(initial, field_turbulence=0.0, field_exchange=0.0)
    active = [
        index
        for index, amplitude in enumerate(initial.sparse_amplitude)
        if amplitude > 0.0
    ]
    if not active:
        return initial
    targets = []
    for cell in cells:
        modal, _ = reference_modal_residual(cell, projection_model)
        reference = fixed_rate(cell.reference)
        onset = detect_impact_onset(reference.samples, reference.sample_rate)
        targets.append(modal_power_features(modal, onset))

    columns = []
    for mode in active:
        amplitudes = [0.0] * len(initial.sparse_amplitude)
        amplitudes[mode] = 1.0
        basis_fit = replace(
            projection_model,
            sparse_amplitude=tuple(amplitudes),
            direct_gain=0.0,
            field_gain=1.0,
            output_gain=1.0,
        )
        rendered = []
        for cell in cells:
            candidate = render_sparse(cell, basis_fit, 2.5)
            onset = detect_impact_onset(candidate.samples, candidate.sample_rate)
            rendered.append(modal_power_features(candidate, onset))
        columns.append(np.concatenate(rendered))

    matrix = np.column_stack(columns)
    target = np.concatenate(targets)
    # Relative-energy weighting prevents the loudest band alone from owning
    # the solution, while a 60 dB power floor avoids fitting inaudible bins.
    floor = max(float(np.max(target)) * 1.0e-6, np.finfo(float).tiny)
    row_weight = 1.0 / np.sqrt(np.maximum(target, floor))
    weighted_matrix = matrix * row_weight[:, None]
    weighted_target = target * row_weight
    column_scale = np.maximum(
        np.linalg.norm(weighted_matrix, axis=0), np.finfo(float).tiny
    )
    normalized = weighted_matrix / column_scale
    normalized_coefficients = nnls(normalized, weighted_target)[0]
    power_coefficients = normalized_coefficients / column_scale

    amplitudes = [0.0] * len(initial.sparse_amplitude)
    for active_index, mode in enumerate(active):
        power = max(power_coefficients[active_index], 0.0)
        amplitudes[mode] = float(np.sqrt(power))
    gain = float(np.linalg.norm(amplitudes))
    if gain > 0.0:
        amplitudes = [amplitude / gain for amplitude in amplitudes]
    return replace(
        initial,
        sparse_amplitude=tuple(amplitudes),
        field_gain=float(np.clip(gain, 0.0, 4.0)),
    )


def reference_modal_residual(
    cell: CrashFitCell, fit: CrashFit
) -> tuple[AudioBuffer, AudioBuffer]:
    reference = fixed_rate(cell.reference)
    onset = detect_impact_onset(reference.samples, reference.sample_rate)
    start = min(reference.samples.size, onset + round(0.2 * reference.sample_rate))
    end = min(reference.samples.size, start + round(4.0 * reference.sample_rate))
    modes = tuple(
        DampedMode(frequency, decay / np.log(1000.0), 1.0, 0.0)
        for frequency, amplitude in zip(
            fit.sparse_frequency_hz,
            fit.sparse_amplitude,
        )
        for decay in (_body_t60(fit, frequency, reference.sample_rate),)
        if amplitude > 0.0
    )
    modal_samples = np.zeros_like(reference.samples)
    if modes and end > start:
        fitted = refit_mode_amplitudes(
            reference.samples[start:end], reference.sample_rate, modes
        )
        modal_samples[start:] = resynthesize_modes(
            fitted, reference.sample_rate, reference.samples.size - start
        )
        fade_samples = min(
            round(0.02 * reference.sample_rate), modal_samples.size - start
        )
        if fade_samples > 0:
            phase = np.linspace(0.0, 0.5 * np.pi, fade_samples, endpoint=True)
            modal_samples[start : start + fade_samples] *= np.square(np.sin(phase))
    return (
        AudioBuffer(modal_samples, reference.sample_rate),
        AudioBuffer(reference.samples - modal_samples, reference.sample_rate),
    )


def _body_t60(fit: CrashFit, frequency_hz: float, sample_rate: int) -> float:
    """Match the production ERB-frequency/log-T60 interpolation."""
    return float(
        interpolate_t60(
            frequency_hz,
            (0.0, *fit.body_decay_frequency_hz, 0.5 * sample_rate),
            fit.body_decay_seconds,
            (True, *fit.body_decay_active, True),
            sample_rate,
        )
    )


def modal_passes() -> tuple[EspritPass, ...]:
    return tuple(
        EspritPass(minimum, maximum, count, 384, 6144, 257)
        for minimum, maximum, count in (
            (100.0, 700.0, 8),
            (600.0, 1400.0, 10),
            (1300.0, 3000.0, 12),
            (2800.0, 6000.0, 14),
            (5500.0, 9000.0, 16),
            (8500.0, 12000.0, 18),
            (11500.0, 14000.0, 18),
            (13500.0, 15000.0, 14),
        )
    )


def stratified_modes(
    modes: tuple[DampedMode, ...], count: int
) -> tuple[DampedMode, ...]:
    """Keep energetic ridges without allowing the low spectrum to use every slot."""
    if count < 1:
        return ()
    score = lambda mode: mode.amplitude * np.sqrt(mode.decay_seconds)
    selected: list[DampedMode] = []
    # Equal quotas keep paintable packet centres across the audible cymbal
    # spectrum. Remaining slots are filled globally when a band is sparse.
    bands = (
        (100.0, 700.0, 4),
        (700.0, 1500.0, 4),
        (1500.0, 3500.0, 4),
        (3500.0, 7000.0, 4),
        (7000.0, 10000.0, 4),
        (10000.0, 12500.0, 2),
        (12500.0, 15000.1, 2),
    )
    scale = count / sum(quota for _, _, quota in bands)
    quotas = [max(0, round(quota * scale)) for _, _, quota in bands]
    while sum(quotas) < count:
        quotas[sum(quotas) % len(quotas)] += 1
    while sum(quotas) > count:
        index = max(range(len(quotas)), key=quotas.__getitem__)
        quotas[index] -= 1
    for (minimum, maximum, _), quota in zip(bands, quotas):
        candidates = sorted(
            (mode for mode in modes if minimum <= mode.frequency_hz < maximum),
            key=score,
            reverse=True,
        )
        selected.extend(candidates[:quota])
    if len(selected) < count:
        selected_ids = {id(mode) for mode in selected}
        remaining = sorted(
            (mode for mode in modes if id(mode) not in selected_ids),
            key=score,
            reverse=True,
        )
        selected.extend(remaining[: count - len(selected)])
    return tuple(selected[:count])


def stable_modes(
    hits: tuple[np.ndarray, ...],
    sample_rate: int,
    passes: tuple[EspritPass, ...],
    required_hits: int,
) -> tuple[DampedMode, ...]:
    evidence = estimate_modal_evidence(hits, sample_rate, passes, merge_cents=18.0)
    return tuple(
        item.mode
        for item in evidence
        if item.hit_count >= required_hits
        and item.observation_count >= required_hits
        and item.frequency_std_cents <= 12.0
        and item.decay_std_log <= 1.0
    )


def modal_tail(audio: AudioBuffer) -> np.ndarray:
    mono = fixed_rate(audio)
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    start = min(mono.samples.size, onset + round(0.2 * mono.sample_rate))
    return mono.samples[start:]
