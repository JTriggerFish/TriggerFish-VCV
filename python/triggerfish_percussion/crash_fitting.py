"""Coarse perceptual features and staged crash fitting helpers."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import differential_evolution, least_squares

from .alignment import detect_impact_onset
from .audio_io import AudioBuffer, resample_audio
from .crash_model import CrashFit, render_crash
from .descriptors import spectral_trajectories
from .erb import ErbFilterbank
from .transforms import StftConfig, stft


@dataclass(frozen=True)
class CrashFitCell:
    label: str
    reference: AudioBuffer
    strength: float
    location: float = 1.0
    hardness: float = 0.65
    seed: int = 1


@dataclass(frozen=True)
class FeatureSet:
    values: np.ndarray
    names: tuple[str, ...]


REGIONS = (
    ("contact", 0.000, 0.015),
    ("bloom", 0.015, 0.250),
    ("early", 0.250, 1.500),
    ("tail-a", 1.500, 4.000),
    ("tail-b", 4.000, 8.000),
)


def perceptual_features(audio: AudioBuffer) -> FeatureSet:
    mono = audio.mono()
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    signal = mono.samples[onset:]
    transform = stft(signal, mono.sample_rate, StftConfig(2048, 512, 4096))
    bank = ErbFilterbank.create(transform.frequencies_hz, 80.0, 18000.0, 24)
    band_power = bank.apply_power(transform.power)
    density_bins = (transform.frequencies_hz >= 300.0) & (
        transform.frequencies_hz <= 18000.0
    )
    trajectories = spectral_trajectories(
        transform.frequencies_hz[density_bins], transform.power[density_bins]
    )
    values: list[float] = []
    names: list[str] = []
    floor = max(float(np.max(band_power)) * 1.0e-10, np.finfo(float).tiny)
    for region, start, end in REGIONS:
        selected = (transform.times_seconds >= start) & (transform.times_seconds < end)
        if not np.any(selected):
            region_power = np.full(band_power.shape[0], floor)
        else:
            region_power = np.mean(band_power[:, selected], axis=1)
        total = max(float(np.sum(region_power)), floor)
        level_db = 10.0 * np.log10(total)
        shape_db = 10.0 * np.log10(np.maximum(region_power / total, 1.0e-10))
        values.append(level_db / 6.0)
        names.append(f"level/{region}")
        values.extend((shape_db / 12.0).tolist())
        names.extend(f"erb/{region}/{frequency:.0f}" for frequency in bank.centers_hz)
        if np.any(selected):
            density_values = (
                np.log(max(float(np.median(trajectories.centroid_hz[selected])), 1.0))
                / 0.5,
                np.log(max(float(np.median(trajectories.flatness[selected])), 1.0e-8)),
                np.log(max(float(np.median(trajectories.crest[selected])), 1.0e-8)),
                float(np.median(trajectories.flux[selected])) / 0.05,
            )
        else:
            density_values = (0.0, 0.0, 0.0, 0.0)
        values.extend(density_values)
        names.extend(
            f"{descriptor}/{region}"
            for descriptor in ("centroid", "flatness", "crest", "flux")
        )
    return FeatureSet(np.asarray(values), tuple(names))


def fit_level_curve(cells: tuple[CrashFitCell, ...], initial: CrashFit) -> CrashFit:
    strengths = np.asarray([cell.strength for cell in cells])
    target = np.asarray([_early_rms(cell.reference) for cell in cells])
    model = np.asarray(
        [_early_rms(_render_cell(cell, initial, 0.25)) for cell in cells]
    )
    selected = (strengths > 0) & (target > 0) & (model > 0)
    x = np.log(strengths[selected])
    target_slope, target_offset = np.polyfit(x, np.log(target[selected]), 1)
    model_slope, model_offset = np.polyfit(x, np.log(model[selected]), 1)
    gamma = np.clip(initial.strength_gamma + target_slope - model_slope, 0.25, 4.0)
    gain = initial.output_gain * np.exp(target_offset - model_offset)
    return replace(initial, strength_gamma=float(gamma), output_gain=float(gain))


def fit_body(
    cells: tuple[CrashFitCell, ...],
    initial: CrashFit,
    maximum_evaluations: int = 800,
    seed: int = 1234,
) -> tuple[CrashFit, dict[str, float]]:
    names = (
        "resonance_tune",
        "low_decay_scale",
        "middle_decay_scale",
        "high_decay_scale",
        "resonator_coupling",
        "resonator_shift_scale",
        "dispersion_feedback",
        "dispersion_drive",
        "dispersion_excursion_samples",
        "dispersion_low_decay_seconds",
        "dispersion_middle_decay_seconds",
        "dispersion_high_decay_seconds",
        "colour_frequency_hz",
        "colour_gain_db",
        "high_cut_hz",
    )
    bounds = (
        (0.65, 1.65),
        (0.35, 2.2),
        (0.35, 2.2),
        (0.25, 2.2),
        (0.0, 1.0),
        (0.0, 4.0),
        (0.70, 0.995),
        (0.0, 8.0),
        (0.0, 6.0),
        (0.15, 3.0),
        (0.12, 2.0),
        (0.08, 1.2),
        (1800.0, 11000.0),
        (-12.0, 12.0),
        (7000.0, 22000.0),
    )
    targets = tuple(perceptual_features(_fixed_rate(cell.reference)) for cell in cells)

    def residual(vector: np.ndarray) -> np.ndarray:
        fit = replace(initial, **dict(zip(names, vector.tolist())))
        errors = []
        for cell, target in zip(cells, targets):
            rendered = _render_cell(cell, fit, 8.2)
            errors.append(perceptual_features(rendered).values - target.values)
        return np.concatenate(errors)

    def scalar(vector: np.ndarray) -> float:
        error = residual(vector)
        return float(np.mean(np.square(error)))

    population = max(4, maximum_evaluations // (12 * len(names)))
    iterations = max(1, maximum_evaluations // (population * len(names)) - 1)
    generation = 0

    def progress(_: np.ndarray, convergence: float) -> bool:
        nonlocal generation
        generation += 1
        print(
            f"crash global fit generation {generation}/{iterations}, "
            f"convergence {convergence:.4g}",
            flush=True,
        )
        return False

    global_fit = differential_evolution(
        scalar,
        bounds,
        seed=seed,
        popsize=population,
        maxiter=iterations,
        polish=False,
        updating="immediate",
        workers=1,
        callback=progress,
    )
    remaining = max(maximum_evaluations - global_fit.nfev, 24)
    polished = least_squares(
        residual,
        global_fit.x,
        bounds=np.asarray(bounds).T,
        max_nfev=remaining,
        verbose=0,
    )
    fit = replace(initial, **dict(zip(names, polished.x.tolist())))
    diagnostics = {
        "global_objective": float(global_fit.fun),
        "final_objective": float(np.mean(np.square(polished.fun))),
        "global_evaluations": int(global_fit.nfev),
        "polish_evaluations": int(polished.nfev),
    }
    return fit, diagnostics


def _fixed_rate(audio: AudioBuffer) -> AudioBuffer:
    return resample_audio(audio.mono(), 48000)


def _render_cell(cell: CrashFitCell, fit: CrashFit, seconds: float) -> AudioBuffer:
    return AudioBuffer(
        render_crash(
            fit,
            seconds,
            48000,
            cell.strength,
            cell.location,
            cell.hardness,
            cell.seed,
        ),
        48000,
    )


def _early_rms(audio: AudioBuffer) -> float:
    mono = _fixed_rate(audio)
    onset = detect_impact_onset(mono.samples, mono.sample_rate)
    end = min(mono.samples.size, onset + round(0.1 * mono.sample_rate))
    return float(np.sqrt(np.mean(np.square(mono.samples[onset:end]))))
