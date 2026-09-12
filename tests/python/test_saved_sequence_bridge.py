"""Opt-in integration with the existing developer workbench server and Wasm."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

pytestmark = pytest.mark.skipif(
    os.environ.get("TF_TEST_WORKBENCH_BRIDGE") != "1",
    reason="Requires built workbench, running reference server and EMSDK_NODE",
)
ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def voice():
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", ROOT)
    fit = json.loads(
        (ROOT / "workbench/web/crash_calibration.fit.json").read_text(encoding="utf8")
    )
    fit["controls"]["event"].update(strength=0.31, location=0.25, seed=2917)
    # Override both the loaded patch's routing and its output gain.
    for connection in fit["instrument"]["connections"]:
        if connection["from"] == "body.audio":
            connection["enabled"] = False
    fit["instrument"]["nodes"][-1]["parameters"]["model_level_db"] = -18
    try:
        yield renderer, fit
    finally:
        renderer.close()


@pytest.mark.parametrize("override", [{}, {"strength": 0.47, "seed": 271}])
def test_saved_sequence_matches_snapshot_with_edited_routing_and_gesture(
    voice, override
):
    renderer, fit = voice
    saved = SavedFitRenderer(renderer, fit)
    expected = saved.render(saved.initial, 0.5, event=override)
    actual = renderer.decode(
        renderer.request(
            command="renderSequence",
            fit=fit,
            seconds=0.5,
            hits=[dict(time=0, **override)],
        )["pcm"]
    )
    assert np.max(np.abs(expected)) > 0
    np.testing.assert_array_equal(actual, expected)


def test_saved_sequence_rejects_ambiguous_configuration(voice):
    renderer, fit = voice
    with pytest.raises(RuntimeError, match="not both"):
        renderer.request(command="renderSequence", fit=fit, parameters={}, hits=[])
