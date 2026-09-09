"""Broad crash fitting edits cannot become independent ridge/gain fitting."""

import sys
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip("torch")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from refine_crash_ring_balance import broad_start, KNOTS
from triggerfish_percussion.coarse_observation_fit import interpolation_weights


def test_middle_warp_is_ordered_and_leaves_low_high_and_dynamics_alone():
    base = dict(
        model_level_db=-1,
        body_excitation=4,
        body_decay_seconds_0=20,
        field_beat_depth=0.2,
        field_beat_rate_tilt=0.25,
    )
    frequencies = np.geomspace(125, 15000, 24)
    for i in range(32):
        base[f"resolved_frequency_{i}"] = float(frequencies[min(i, 23)])
        base[f"resolved_level_{i}"] = -12 if i < 24 else -72
    for scale in (0.9, 1, 1.1):
        p, amplitudes = broad_start(base, scale)
        f = [p[f"resolved_frequency_{i}"] for i in range(24)]
        assert np.all(np.diff(f) > 0)
        assert f[0] == frequencies[0] and f[-1] == frequencies[-1]
        actual = 10 ** (np.array([p[f"resolved_level_{i}"] for i in range(24)]) / 20)
        np.testing.assert_allclose(actual, interpolation_weights(f, KNOTS) @ amplitudes)
        for key, value in base.items():
            if not key.startswith(("resolved_level_", "resolved_frequency_")):
                assert p[key] == value
        assert p["resolved_level_31"] == -72
