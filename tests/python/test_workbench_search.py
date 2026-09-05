"""Inactive controls must remain fixed, rather than drifting to bounds."""

from types import SimpleNamespace

import numpy as np

from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.workbench_multistart import refine_candidate_starts


def test_raw_bad_start_can_win_after_short_refit(tmp_path, monkeypatch):
    renderer = SimpleNamespace(
        initial={"gain": 0},
        metadata={"descriptors": [dict(key="gain", minimum=0, maximum=1)]},
        render=lambda parameters, seconds: np.array([parameters["gain"]]),
    )
    loss = SimpleNamespace(
        residual=lambda audio, regions=range(5): 10 * (audio - 0.4),
        diagnostics=lambda audio: {},
    )
    search = Search(renderer, loss, tmp_path)
    monkeypatch.setattr(Search, "save", lambda self: None)
    records = refine_candidate_starts(
        search,
        [("poor raw gain", {"gain": 1})],
        lambda values: {"gain": (0, 1)},
        count=1,
    )
    assert abs(search.parameters["gain"] - 0.4) < 0.001
    assert records[0]["raw_score"] > search.history[0]["before"]
    assert records[0]["refined_score"] < search.history[0]["before"]


def test_kick_fit_uses_modes_not_eq_or_individual_damping():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "tools/kick_fit_stages.py"
    spec = importlib.util.spec_from_file_location("kick_fit_stages", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    parameters = {f"resonance_level_{i}": -72 for i in range(16)}
    parameters.update(resonance_level_0=0, resonance_frequency_0=55)
    bounds = dict(
        item for _, stage in module.stages_for(parameters) for item in stage.items()
    )
    assert "resonance_frequency_0" in bounds
    assert "resonance_level_0" not in bounds  # fixed reference prominence
    assert "resonance_frequency_1" not in bounds
    assert not any(
        key.startswith(("band_", "equalizer", "colour_", "low_cut", "high_cut"))
        for key in bounds
    )
    assert not any(
        key.startswith(("resonance_centre_", "resonance_edge_")) for key in bounds
    )
    assert {key for key in bounds if key.startswith("resonance_decay")} == {
        "resonance_decay_seconds",
        "resonance_decay_tilt",
    }
    parameters.update(resonance_level_1=-6, resonance_frequency_1=1200)
    bounds = dict(
        item for _, stage in module.stages_for(parameters) for item in stage.items()
    )
    assert "resonance_level_0" not in bounds
    assert "resonance_level_1" in bounds
    assert bounds["resonance_frequency_1"][1] >= 1200


def test_fit_recovers_active_control_and_freezes_dead_direction(tmp_path, monkeypatch):
    renderer = SimpleNamespace(
        initial={"active": 0.0, "dead": 0.7},
        metadata={
            "descriptors": [
                {"key": key, "minimum": -1, "maximum": 1} for key in ("active", "dead")
            ]
        },
        render=lambda parameters, seconds: np.array([parameters["active"]]),
    )
    loss = SimpleNamespace(
        residual=lambda audio, regions=range(5): 10 * (audio - 0.4),
        diagnostics=lambda audio: {},
    )
    search = Search(renderer, loss, tmp_path)
    monkeypatch.setattr(search, "save", lambda: None)
    search.stage("fixture", {"active": (-100, 100), "dead": (-100, 100)}, 10)
    assert abs(search.parameters["active"] - 0.4) < 0.001
    assert search.parameters["dead"] == 0.7


def test_narrow_bounds_cannot_accept_regression_from_actual_patch(
    tmp_path, monkeypatch
):
    renderer = SimpleNamespace(
        initial={"active": 0.9},
        metadata={"descriptors": [{"key": "active", "minimum": 0, "maximum": 1}]},
        render=lambda parameters, seconds: np.array([parameters["active"]]),
    )
    loss = SimpleNamespace(
        residual=lambda audio, regions: 10 * (audio - 0.9), diagnostics=lambda audio: {}
    )
    search = Search(renderer, loss, tmp_path)
    monkeypatch.setattr(search, "save", lambda: None)
    search.stage("narrower box", {"active": (0, 0.5)}, 10)
    assert search.parameters["active"] == 0.9
    assert not search.history[-1]["selected"]
