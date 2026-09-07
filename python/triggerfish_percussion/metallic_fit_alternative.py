"""Refit a changed state regime before comparing it to the current sound."""

import numpy as np

from .observation_fit_polish import polish_observation
from .workbench_search import Search
from .modal_fit_initialization import spectral_mode_candidates
from .fit_rerender import load_checked_start


def late_ridge_parameters(parameters, reference, rate, low=250):
    """Stable tail peaks as coherent handles, plus three diffuse upper handles."""
    peaks = spectral_mode_candidates(
        reference, rate, count=28, low=low, high=3500, start=0.4, end=3
    )
    if not peaks:
        raise ValueError("No resolved tail peaks for a coherent-body alternative")
    values = dict(parameters, bloom_rate=0, field_exchange=0, bloom_phase_diffusion=0)
    maximum = max(p["power_db"] for p in peaks)
    for i in range(32):
        values[f"resolved_level_{i}"] = -72
        if i < len(peaks):
            values[f"resolved_frequency_{i}"] = peaks[i]["frequency"]
            values[f"resolved_level_{i}"] = max(-42, peaks[i]["power_db"] - maximum)
            values[f"resolved_turbulence_{i}"] = 0
    for i, frequency in enumerate((6000.0, 9000.0, 14000.0), len(peaks)):
        values[f"resolved_frequency_{i}"] = frequency
        values[f"resolved_level_{i}"] = -18
        values[f"resolved_turbulence_{i}"] = 1
    return values


def refit_alternative(search, values, name, iterations=8):
    before = float(np.linalg.norm(search.residual(search.parameters)))
    directory = search.output / name
    directory.mkdir(parents=True, exist_ok=True)
    trial = Search(
        search.renderer,
        search.loss,
        directory,
        search.seconds,
        search.name,
        search.seeds,
    )
    trial.parameters = dict(values)
    if (directory / "search.json").exists():
        saved, provenance = load_checked_start(
            search.renderer, directory / "search.json"
        )
        trial.parameters, trial.history = saved["parameters"], saved["history"]
        trial.history.append(
            dict(stage="checked alternative warm start", provenance=provenance)
        )
    bounds = dict(
        field_gain=(0.1, 4),
        body_decay_seconds_0=(0.1, 30),
        body_decay_seconds_7=(0.02, 15),
    )
    if trial.parameters.get("body_decay_active_1", 0) > 0.5:
        bounds["body_decay_seconds_1"] = (5, 30)
    trial.stage("alternative level and damping", bounds, 8)
    polish_observation(trial, iterations)
    after = float(np.linalg.norm(search.residual(trial.parameters)))
    selected = after < before
    if selected:
        search.parameters = dict(trial.parameters)
    search.history.append(
        dict(
            stage=name,
            before=before,
            after=after,
            selected=selected,
            trial_parameters=trial.parameters,
            parameters=dict(search.parameters),
        )
    )
    search.save()
