"""The developer scalar fitter keeps valid, best-evaluated visible controls."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "metal_refinement", Path(__file__).parents[2] / "tools/refine_metal_perceptual.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fake_search(tmp_path, initial):
    class Loss:
        specification = {"implementation": "known quadratic"}

        def score(self, samples):
            return float((samples[0] - 0.7) ** 2)

    return SimpleNamespace(
        renderer=SimpleNamespace(
            metadata={"descriptors": [{"key": "position", "minimum": 0, "maximum": 1}]}
        ),
        parameters={"position": initial, "fixed": 42},
        loss=Loss(),
        output=tmp_path,
        history=[],
        audio=lambda p: np.array([p["position"]]),
        save=lambda: None,
    )


def test_scalar_fit_recovers_target_and_preserves_fixed_values(tmp_path):
    search = fake_search(tmp_path, 0.1)
    module.scalar_stage(search, "known target", {"position": (0, 1)}, 50)
    assert search.parameters["position"] == pytest.approx(0.7, abs=0.005)
    assert search.parameters["fixed"] == 42
    assert search.history[0]["after"] < search.history[0]["before"]


def test_narrow_box_cannot_authorize_regression(tmp_path):
    search = fake_search(tmp_path, 0.7)
    module.scalar_stage(search, "exclude optimum", {"position": (0, 0.2)}, 30)
    assert search.parameters["position"] == 0.7
    assert search.history[0]["after"] == search.history[0]["before"]


def test_never_search_outside_ui_range(tmp_path):
    search = fake_search(tmp_path, 0.1)
    with pytest.raises(ValueError, match="exposed control range"):
        module.scalar_stage(search, "invalid", {"position": (-1, 1)}, 30)
