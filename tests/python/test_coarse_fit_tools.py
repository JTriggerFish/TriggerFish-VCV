"""Sparse-curve validation and rejection of a regressive joint-fit result."""

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("torch")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
import fit_joint_coarse_decay as joint
from refine_coarse_metal import damping_variant


@pytest.mark.parametrize(
    "frequencies",
    [[40], [15000], [np.nan], [800, 400], [500, 500], list(range(100, 800, 100))],
)
def test_sparse_curve_rejects_invalid_knots(frequencies):
    with pytest.raises(ValueError):
        damping_variant({}, frequencies)


def test_joint_fit_preserves_original_when_new_curve_is_worse(monkeypatch, tmp_path):
    original = {"body_decay_seconds_0": 4.0, "body_decay_seconds_7": 1.0}
    for i in range(32):
        original[f"resolved_level_{i}"] = -20.0 if i < 6 else -72.0
        original[f"resolved_frequency_{i}"] = 100.0 * (i + 1)
    # An existing interior knot is intentionally not on the endpoint curve.
    original.update(
        body_decay_active_1=1, body_decay_frequency_1=400, body_decay_seconds_1=10.0
    )
    search = SimpleNamespace(
        parameters=original.copy(),
        output=tmp_path,
        history=[],
        save=lambda: None,
        residual=lambda p: np.array([p["body_decay_seconds_1"] - 10.0]),
    )
    monkeypatch.setattr(
        joint,
        "least_squares",
        lambda fun, start, **kwargs: SimpleNamespace(
            x=start, message="test result", nfev=1
        ),
    )
    joint.fit(search, [400], 1)
    assert search.parameters == original
    assert not search.history[-1]["accepted"]
    assert search.history[-1]["after"] == 0
