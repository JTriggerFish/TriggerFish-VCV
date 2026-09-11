"""Refinement must reproduce the saved gesture, not the factory target gesture."""

from copy import deepcopy
from types import SimpleNamespace

import pytest

from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer


def fixture():
    reference = dict(sha256="source", referenceGainDb=-6, sampleRate=44100)
    fit = dict(
        id="user-edit",
        name="Original",
        reference=reference,
        instrument=dict(
            recipe="metal.cymbal.v1",
            nodes=[dict(parameters=dict(bloom_rate=4))],
            edges=[dict(gain=0.7)],
        ),
        controls=dict(event=dict(strength=0.93, seed=1695, implement=0.5)),
    )
    r = SimpleNamespace(
        sample_rate=44100,
        initial=dict(bloom_rate=2),
        metadata=dict(reference=reference, recipeKey="metal.cymbal.v1"),
        decode=lambda x: x,
    )
    r.request = lambda **kwargs: dict(pcm=kwargs)
    return r, fit


def test_saved_gesture_and_routing_survive_parameter_edit():
    r, fit = fixture()
    original = deepcopy(fit)
    saved = SavedFitRenderer(r, fit)
    result = saved.render(dict(bloom_rate=3), 6)
    actual = result["fit"]
    assert actual["controls"] == fit["controls"]
    assert actual["instrument"]["edges"] == fit["instrument"]["edges"]
    assert actual["instrument"]["nodes"][0]["parameters"] == dict(bloom_rate=3)
    assert fit == original


def test_holdout_seed_does_not_change_the_saved_strike():
    r, fit = fixture()
    saved = SavedFitRenderer(r, fit)
    event = saved.render(saved.initial, 6, seed=42)["fit"]["controls"]["event"]
    assert event == dict(strength=0.93, seed=42, implement=0.5)
    assert fit["controls"]["event"]["seed"] == 1695


def test_snapshot_tracks_parent_without_overwriting_source():
    r, fit = fixture()
    candidate = SavedFitRenderer(r, fit).snapshot(dict(bloom_rate=3), "Candidate")
    assert candidate["id"] != fit["id"]
    assert candidate["parentId"] == fit["id"]
    assert fit["name"] == "Original"


@pytest.mark.parametrize("key,value", [("sha256", "different"), ("referenceGainDb", 0)])
def test_wrong_reference_rejected(key, value):
    r, fit = fixture()
    fit = deepcopy(fit)
    fit["reference"][key] = value
    with pytest.raises(ValueError, match="reference"):
        SavedFitRenderer(r, fit)


def test_wrong_reference_onset_rejected():
    r, fit = fixture()
    fit = deepcopy(fit)
    fit["reference"]["cell"] = dict(onset_seconds=0.1)
    with pytest.raises(ValueError, match="onset"):
        SavedFitRenderer(r, fit)


def test_reference_without_corpus_cell_is_valid():
    r, fit = fixture()
    fit = deepcopy(fit)
    fit["reference"]["cell"] = None
    assert SavedFitRenderer(r, fit).initial == dict(bloom_rate=4)
