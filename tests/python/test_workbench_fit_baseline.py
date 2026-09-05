"""A model rebuild must not silently replace the original audition baseline."""

from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from triggerfish_percussion.workbench_fit_baseline import original_baseline


def renderer():
    return SimpleNamespace(
        sample_rate=200,
        initial={},
        metadata=dict(
            reference=dict(sha256="fixture", sampleRate=200, referenceGainDb=2),
            event=dict(strength=0.5),
            rendererSha256="first-model",
        ),
        render=lambda parameters, seconds: np.array([0.1, 0.2]),
    )


def test_model_revision_keeps_baseline(tmp_path):
    engine = renderer()
    first = original_baseline(engine, tmp_path, 0.01)
    engine.metadata["rendererSha256"] = "second-model"
    engine.render = lambda parameters, seconds: np.ones(2)
    assert np.allclose(original_baseline(engine, tmp_path, 0.01), first)


def test_source_gain_and_duration_changes_are_rejected(tmp_path):
    engine = renderer()
    original_baseline(engine, tmp_path, 0.01)
    engine.metadata = deepcopy(engine.metadata)
    engine.metadata["reference"]["referenceGainDb"] = 6
    with pytest.raises(ValueError, match="referenceGainDb"):
        original_baseline(engine, tmp_path, 0.01)
    engine.metadata["reference"]["referenceGainDb"] = 2
    with pytest.raises(ValueError, match="duration"):
        original_baseline(engine, tmp_path, 0.02)


def test_onset_changes_are_rejected_even_for_the_same_source(tmp_path):
    engine = renderer()
    original_baseline(engine, tmp_path, 0.01)
    engine.metadata["reference"]["cell"] = {"onset_seconds": 0.002}
    with pytest.raises(ValueError, match="onset"):
        original_baseline(engine, tmp_path, 0.01)
