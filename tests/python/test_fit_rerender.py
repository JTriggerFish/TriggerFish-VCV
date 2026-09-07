"""Engine changes must not silently retarget an archived fitting experiment."""

import json
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest
from triggerfish_percussion.fit_rerender import (
    load_checked_start,
    prepare_current_render,
)
from triggerfish_percussion.kick_quality_checks import heldout_seeds


def fixture(tmp_path):
    metadata = dict(
        recipeKey="drum.kick.v1",
        rendererSha256="old",
        descriptors=[dict(key="gain")],
        event=dict(seed=1449),
        reference=dict(
            sha256="source",
            sampleRate=44100,
            referenceGainDb=2,
            cell=dict(onset_seconds=0.001),
        ),
    )
    saved = dict(
        metadata=deepcopy(metadata),
        parameters=dict(gain=0.5),
        duration_seconds=1,
        training_seeds=[1449, 1450],
    )
    (tmp_path / "search.json").write_text(json.dumps(saved))
    metadata["rendererSha256"] = "new"
    calls = []
    renderer = SimpleNamespace(
        metadata=metadata,
        initial=dict(gain=0.9),
        sample_rate=44100,
        render=lambda *args: calls.append(args) or np.zeros(44100),
        request=lambda **kwargs: dict(fit={"parameters": kwargs["parameters"]}),
    )
    return renderer, saved, calls


@pytest.mark.parametrize(
    "field",
    [
        "sha256",
        "sampleRate",
        "referenceGainDb",
        "onset",
        "event",
        "recipe",
        "descriptors",
        "parameters",
    ],
)
def test_retargeting_rejected_before_render_or_write(tmp_path, field):
    renderer, saved, calls = fixture(tmp_path)
    if field == "onset":
        renderer.metadata["reference"]["cell"]["onset_seconds"] = 0.1
    elif field == "event":
        renderer.metadata["event"]["seed"] = 5
    elif field == "recipe":
        renderer.metadata["recipeKey"] = "different"
    elif field == "descriptors":
        renderer.metadata["descriptors"] = []
    elif field == "parameters":
        renderer.initial = {"different": 1}
    else:
        renderer.metadata["reference"][field] = "changed"
    output = tmp_path / "output"
    with pytest.raises(ValueError):
        prepare_current_render(renderer, tmp_path, output)
    assert not calls and not output.exists()


def test_rerender_preserves_origin_and_frozen_parameters(tmp_path):
    renderer, saved, calls = fixture(tmp_path)
    output = tmp_path / "output"
    prepare_current_render(renderer, tmp_path, output)
    result = json.loads((output / "search.json").read_text())
    assert result["rerendered_from"]["metadata"] == saved["metadata"]
    assert result["metadata"]["rendererSha256"] == "new"
    assert calls[0][0] == {"gain": 0.5}
    # Changing today's preset cannot change the archived audit baseline.
    first, provenance = load_checked_start(renderer, tmp_path / "search.json")
    renderer.initial["gain"] = 0.1
    second, again = load_checked_start(renderer, tmp_path / "search.json")
    assert first == second == saved and provenance == again


def test_heldout_excludes_every_training_seed_and_wraps():
    saved = dict(
        metadata=dict(event=dict(seed=1449)), training_seeds=[1449, 1450, 1452]
    )
    assert heldout_seeds(saved) == [1451, 1453, 1454]
    saved.update(metadata=dict(event=dict(seed=0xFFFFFFFF)), training_seeds=[None, 0])
    assert heldout_seeds(saved) == [1, 2, 3]
    del saved["training_seeds"]
    with pytest.raises(ValueError, match="Training seeds"):
        heldout_seeds(saved)
