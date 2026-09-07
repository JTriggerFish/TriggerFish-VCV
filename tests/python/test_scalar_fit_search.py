"""Known scalar optimum and strict budgets without expensive synthesis."""

from types import SimpleNamespace
import numpy as np
import pytest
from triggerfish_percussion.scalar_fit_search import refine_scalar


def search():
    obj = SimpleNamespace(
        parameters={"gain": 0.1, "pitch_hz": 40},
        history=[],
        evaluations=0,
        seeds=(1, 2),
        name="test",
    )
    obj.renderer = SimpleNamespace(
        metadata={
            "descriptors": [
                dict(key="gain", minimum=0, maximum=1, scale="linear"),
                dict(key="pitch_hz", minimum=20, maximum=200, scale="logarithmic"),
            ]
        }
    )
    obj.loss = SimpleNamespace(
        score=lambda x: (x[0] - 0.7) ** 2 + np.log(x[1] / 80) ** 2
    )

    def audio(values, seed):
        obj.evaluations += 1
        return [values["gain"], values["pitch_hz"]]

    obj.audio = audio
    obj.save = lambda: None
    return obj


def test_known_optimum_and_equal_seed_budget():
    obj = search()
    row = refine_scalar(obj, {"gain": (0, 1), "pitch_hz": (20, 200)}, 100)
    assert row["after"] < 1e-5
    assert row["renders"] == row["parameter_evaluations"] * 2
    assert row["parameter_evaluations"] <= 100
    assert len(row["influence"]) == 2


def test_tiny_budget_retains_valid_candidate():
    obj = search()
    row = refine_scalar(obj, {"gain": (0, 1)}, 2)
    assert row["parameter_evaluations"] == 2
    assert row["after"] <= row["before"]


def test_outside_ui_bounds_rejected():
    with pytest.raises(ValueError, match="Bounds outside UI"):
        refine_scalar(search(), {"gain": (-1, 1)}, 10)


def test_descriptor_scales_allow_zero_hold_and_logarithmic_pitch():
    from triggerfish_percussion.fit_parameter_box import ParameterBox

    box = ParameterBox(
        {"thump_hold_seconds": 0, "pitch_hz": 20},
        {"thump_hold_seconds": (0, 0.08), "pitch_hz": (20, 200)},
        [
            dict(key="thump_hold_seconds", minimum=0, maximum=0.08, scale="linear"),
            dict(key="pitch_hz", minimum=20, maximum=200, scale="logarithmic"),
        ],
    )
    assert box.initial.tolist() == [0, 0]
    middle = box.unpack(np.array([0.5, 0.5]))
    assert middle["thump_hold_seconds"] == pytest.approx(0.04)
    assert middle["pitch_hz"] == pytest.approx(np.sqrt(20 * 200))
