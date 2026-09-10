import numpy as np
import pytest
from triggerfish_percussion.modulation_signature import (
    modulation_signature,
    excess_motion,
)


def test_distinguishes_depth_and_synchrony():
    t = np.arange(6 * 16000) / 16000

    def sound(depth, rates):
        return sum(
            (1 + depth * np.cos(2 * np.pi * r * t)) * np.sin(2 * np.pi * f * t)
            for f, r in zip((130, 240, 400, 750), rates)
        ) * np.exp(-t)

    deep = modulation_signature(sound(0.7, [3] * 4), 16000)
    light = modulation_signature(sound(0.1, [3] * 4), 16000)
    irregular = modulation_signature(sound(0.7, [1.25, 2.5, 4.75, 7.5]), 16000)
    assert deep["bands"][0]["dominant_hz"] == 3
    assert deep["bands"][0]["depth"] > 5 * light["bands"][0]["depth"]
    assert deep["mean_synchrony"] > irregular["mean_synchrony"] + 0.3
    assert deep["common_line_strength"] > 2 * irregular["common_line_strength"]
    assert excess_motion(deep, light) > 0
    assert excess_motion(deep, deep) == 0


def test_slow_and_fast_motion_are_not_interchangeable():
    t = np.arange(6 * 16000) / 16000

    def measure(rate):
        audio = (
            (1 + 0.3 * np.cos(2 * np.pi * rate * t))
            * np.sin(2 * np.pi * 130 * t)
            * np.exp(-t)
        )
        return modulation_signature(audio, 16000)["bands"][0]

    slow, fast = measure(0.75), measure(4)
    assert slow["slow_depth"] > 10 * slow["fast_depth"]
    assert fast["fast_depth"] > 10 * fast["slow_depth"]
    assert slow["dominant_hz"] == 0.75
    assert fast["dominant_hz"] == 4


def test_invalid_and_silent():
    with pytest.raises(ValueError):
        modulation_signature(np.zeros(100), 16000)
    with pytest.raises(ValueError):
        modulation_signature(np.zeros(96000), float("nan"))
    result = modulation_signature(np.zeros(96000), 16000)
    assert result["mean_synchrony"] == 0
    assert all(r["depth"] == 0 for r in result["bands"])
    with pytest.raises(ValueError, match="no measurable"):
        excess_motion(result, result)


def test_flutter_is_not_invisible_to_slow_beating_diagnostic():
    t = np.arange(6 * 16000) / 16000
    audio = (
        (1 + 0.4 * np.cos(2 * np.pi * 30 * t))
        * np.sin(2 * np.pi * 400 * t)
        * np.exp(-t)
    )
    band = modulation_signature(audio, 16000)["bands"][2]
    assert band["flutter_dominant_hz"] == 30
    assert band["flutter_depth"] > 0.2
    assert band["flutter_depth"] > 10 * band["depth"]
