"""Structured crash trials must not mutate geometry, events or shared gains."""

from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("torch")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from refine_crash_beating import CrashObjective, KNOTS, project_levels
from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from screen_stable_crash_texture import allocation_tilt


def test_crash_objective_records_the_actual_weights():
    loss = object.__new__(CrashObjective)
    loss.components = lambda _: dict(mel=2, attack_mel=3, bloom=4, texture=5)
    assert np.isclose(loss.score(None), 2 + 0.3 * 3 + 0.05 * 4 + 0.3 * 5)


def test_project_levels_preserves_all_other_controls_and_smooth_shape():
    source = dict(model_level_db=-12, body_excitation=1.7)
    frequencies = np.geomspace(110, 13500, 32)
    weights = interpolation_weights(frequencies, KNOTS)
    levels = 20 * np.log10(weights @ np.array([0.1, 0.3, 0.2, 0.4, 0.1, 0.15]))
    for i, (frequency, level) in enumerate(zip(frequencies, levels)):
        source[f"resolved_frequency_{i}"] = float(frequency)
        source[f"resolved_level_{i}"] = float(level)
        source[f"resolved_allocation_{i}"] = 1
    target = dict(source, resolved_frequency_12=900, resolved_level_31=-72)
    result = project_levels(target, source)
    for key in target:
        if not key.startswith("resolved_level_"):
            assert result[key] == target[key]
    assert result["resolved_level_31"] == -72
    active = [i for i in range(32) if target[f"resolved_level_{i}"] > -72]
    expected = interpolation_weights(
        [target[f"resolved_frequency_{i}"] for i in active], KNOTS
    )
    actual = 10 ** (np.array([result[f"resolved_level_{i}"] for i in active]) / 20)
    np.testing.assert_allclose(
        actual, expected @ np.array([0.1, 0.3, 0.2, 0.4, 0.1, 0.15]), rtol=1e-10
    )


def test_allocation_tilt_only_changes_active_allocation_weights():
    source = {f"resolved_level_{i}": -12 if i < 3 else -72 for i in range(32)}
    source.update({f"resolved_frequency_{i}": 120 * (i + 1) ** 2 for i in range(32)})
    source.update({f"resolved_allocation_{i}": 1 for i in range(32)})
    source["field_phase_bandwidth"] = 0
    assert allocation_tilt(source, 0) == source
    result = allocation_tilt(source, 1)
    assert (
        0
        < result["resolved_allocation_0"]
        < result["resolved_allocation_1"]
        < result["resolved_allocation_2"]
        == 4
    )
    for key in source:
        if key not in [f"resolved_allocation_{i}" for i in range(3)]:
            assert result[key] == source[key]
