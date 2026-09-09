"""Audition edits remain broad, public, and separate from geometry and gains."""

import numpy as np
import pytest
from pathlib import Path
import sys

pytest.importorskip("torch")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from refine_crash_audition import audition_parameters, BOUNDS


def test_low_shelf_and_shared_edits_preserve_other_controls():
    base = {"model_level_db": -12, "body_excitation": 4, "field_turbulence": 1.5}
    for i in range(32):
        base[f"resolved_frequency_{i}"] = 120 * (i + 1)
        base[f"resolved_level_{i}"] = -12 if i < 24 else -72
    result = audition_parameters(base, np.full(6, 0.5))
    assert result["resolved_level_0"] == -16
    assert result["resolved_level_1"] == -16
    assert -16 < result["resolved_level_2"] < -12
    assert result["resolved_level_3"] == -12
    assert "low_prominence_db" not in result
    for key in base:
        if key not in BOUNDS and not key.startswith("resolved_level_"):
            assert result[key] == base[key]
    assert result["resolved_level_24"] == -72


@pytest.mark.parametrize("position", [np.zeros(5), np.full(6, np.nan), np.full(6, 2)])
def test_invalid_search_coordinates_fail(position):
    with pytest.raises(ValueError):
        audition_parameters({}, position)
