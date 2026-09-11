"""Observation balancing must not alter modal dynamics or harmonic layout."""

import importlib
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("torch")
pytest.importorskip("auraloss")


@pytest.mark.parametrize("absolute", [False, True])
def test_observation_paint_has_no_frequency_specific_lock(monkeypatch, absolute):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[2] / "tools"))
    module = importlib.import_module("refine_gong_layered_timing")
    base = {}
    for i in range(32):
        base[f"resolved_frequency_{i}"] = 120 * (i + 1)
        base[f"resolved_level_{i}"] = -15
    original = dict(base)
    anchors = [360, 900, 2000, 4500, 9000, 14000]
    result = module.paint(base, [-50] * 6, absolute, anchors)
    assert result["resolved_level_0"] < -40
    assert result["resolved_level_1"] < -40
    assert result["resolved_level_2"] < -40
    assert base == original


def test_layered_balance_keeps_centres_and_low_pair(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[2] / "tools"))
    module = importlib.import_module("refine_gong_layered_bloom")
    base = {"body_decay_seconds_0": 10, "bloom_rate": 4, "model_level_db": 0}
    for i in range(32):
        base[f"resolved_frequency_{i}"] = 120 * (i + 1)
        base[f"resolved_level_{i}"] = 6 - i
    base["resolved_level_31"] = -72
    original = dict(base)
    result = module.balanced_parameters(base)
    assert base == original
    assert result.keys() == base.keys()
    assert np.isfinite(list(result.values())).all()
    for key in base:
        if not key.startswith("resolved_level_"):
            assert result[key] == base[key]
    for i in (0, 1):
        assert result[f"resolved_level_{i}"] == base[f"resolved_level_{i}"] - 16
    assert result["resolved_level_31"] == -72
    assert all(-72 <= result[f"resolved_level_{i}"] <= 6 for i in range(32))

    curve = importlib.import_module("polish_gong_layered_observation")
    painted = curve.paint(base, [-16, 0, 20, 40])
    assert painted["resolved_level_31"] == -72
    for key in base:
        if not key.startswith("resolved_level_"):
            assert painted[key] == base[key]
    for i in (0, 1):
        assert painted[f"resolved_level_{i}"] == base[f"resolved_level_{i}"] - 16
