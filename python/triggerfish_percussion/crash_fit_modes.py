"""Persistent-mode extraction, resynthesis, and sparse-bank fitting."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer
from .crash_fit_common import (
    CrashFitCell,
    fixed_rate,
    modal_fit_window,
    render_sparse,
)
from .crash_model import CrashFit
from .modal_evidence import EspritPass, estimate_modal_evidence
from .modes import DampedMode, refit_mode_amplitudes, resynthesize_modes


def fit_sparse_modes(
    cells: tuple[CrashFitCell, ...], initial: CrashFit
) -> tuple[CrashFit, dict[str, object]]:
    passes = modal_passes()
    tails = tuple(modal_tail(cell.reference) for cell in cells)
    required_hits = max(1, int(np.ceil(0.75 * len(cells))))
    stable = stable_modes(tails, 48000, passes, required_hits)
    ranked = sorted(
        stable,
        key=lambda mode: mode.amplitude * np.sqrt(mode.decay_seconds),
        reverse=True,
    )[: len(initial.sparse_frequency_hz)]
    selected = tuple(sorted(ranked, key=lambda mode: mode.frequency_hz))
    count = len(initial.sparse_frequency_hz)
    frequencies = list(initial.sparse_frequency_hz)
    decays = list(initial.sparse_decay_seconds)
    amplitudes = [0.0] * count
    phases = [0.0] * count
    for index, mode in enumerate(selected):
        frequencies[index] = mode.frequency_hz
        decays[index] = np.log(1000.0) * mode.decay_seconds
        amplitudes[index] = mode.amplitude
        phases[index] = mode.phase_radians
    peak = max(amplitudes, default=0.0)
    if peak > 0:
        amplitudes = [amplitude / peak for amplitude in amplitudes]
    fitted = replace(
        initial,
        sparse_frequency_hz=tuple(float(value) for value in frequencies),
        sparse_decay_seconds=tuple(float(value) for value in decays),
        sparse_amplitude=tuple(float(value) for value in amplitudes),
        sparse_phase_radians=tuple(float(value) for value in phases),
        sparse_tune=1.0,
        sparse_decay_scale=1.0,
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
    active = [
        index
        for index, amplitude in enumerate(initial.sparse_amplitude)
        if amplitude > 0.0
    ]
    if not active:
        return initial
    targets = []
    for cell in cells:
        modal, _ = reference_modal_residual(cell, initial)
        reference = fixed_rate(cell.reference)
        onset = detect_impact_onset(reference.samples, reference.sample_rate)
        targets.append(modal_fit_window(modal, onset))

    columns = []
    for mode in active:
        for phase in (0.0, 0.5 * np.pi):
            amplitudes = [0.0] * len(initial.sparse_amplitude)
            phases = [0.0] * len(initial.sparse_phase_radians)
            amplitudes[mode] = 1.0
            phases[mode] = phase
            basis_fit = replace(
                initial,
                sparse_amplitude=tuple(amplitudes),
                sparse_phase_radians=tuple(phases),
                direct_gain=0.0,
                sparse_gain=1.0,
                dense_gain=0.0,
                output_gain=1.0,
            )
            rendered = []
            for cell in cells:
                candidate = render_sparse(cell, basis_fit, 2.5)
                onset = detect_impact_onset(candidate.samples, candidate.sample_rate)
                rendered.append(modal_fit_window(candidate, onset))
            columns.append(np.concatenate(rendered))

    matrix = np.column_stack(columns)
    target = np.concatenate(targets)
    scales = np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).tiny)
    normalized = matrix / scales
    normalized_coefficients = np.linalg.lstsq(normalized, target, rcond=1.0e-5)[0]
    coefficients = normalized_coefficients / scales

    amplitudes = [0.0] * len(initial.sparse_amplitude)
    phases = [0.0] * len(initial.sparse_phase_radians)
    for active_index, mode in enumerate(active):
        cosine = coefficients[2 * active_index]
        sine = coefficients[2 * active_index + 1]
        amplitudes[mode] = float(np.hypot(cosine, sine))
        phases[mode] = float(np.arctan2(sine, cosine))
    gain = float(np.linalg.norm(amplitudes))
    if gain > 0.0:
        amplitudes = [amplitude / gain for amplitude in amplitudes]
    return replace(
        initial,
        sparse_amplitude=tuple(amplitudes),
        sparse_phase_radians=tuple(phases),
        sparse_gain=float(np.clip(gain, 0.0, 4.0)),
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
        for frequency, decay, amplitude in zip(
            fit.sparse_frequency_hz,
            fit.sparse_decay_seconds,
            fit.sparse_amplitude,
        )
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


def modal_passes() -> tuple[EspritPass, ...]:
    return tuple(
        EspritPass(minimum, maximum, count, 384, 6144, 257)
        for minimum, maximum, count in (
            (100.0, 700.0, 8),
            (600.0, 1400.0, 10),
            (1300.0, 3000.0, 12),
            (2800.0, 6000.0, 14),
        )
    )


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
