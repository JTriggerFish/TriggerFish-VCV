"""Fail-closed report publication using the current exact renderer, not file age."""

import json
import numpy as np
from .audio_io import read_wav
from .workbench_fit_baseline import check_reference
from .fit_reference import aligned_reference


def check_saved_snapshot(fit, saved):
    """Silent/unused controls and reference metadata must match, too."""
    metadata = saved["metadata"]
    if fit.get("schema") != "triggerfish.percussion.fit/v1":
        raise ValueError("Unsupported saved UI fit schema")
    check_reference(
        dict(reference=fit["reference"], event=fit["controls"]["event"]), metadata
    )
    instrument = fit["instrument"]
    if instrument["recipe"] != metadata["recipeKey"]:
        raise ValueError("Saved UI fit recipe differs from the search")
    parameters = {}
    for node in instrument["nodes"]:
        for key, value in node["parameters"].items():
            if key in parameters:
                raise ValueError("Duplicate saved UI fit parameter")
            parameters[key] = value
    if parameters != saved["parameters"]:
        raise ValueError("Saved UI fit parameters differ from the search")


def verify_candidate(renderer, directory):
    saved = json.loads((directory / "search.json").read_text(encoding="utf8"))
    metadata = saved["metadata"]
    check_reference(metadata, renderer.metadata)
    for key in ("recipeKey", "rendererSha256", "descriptors"):
        if not metadata.get(key) or metadata[key] != renderer.metadata.get(key):
            raise ValueError(
                f"Stale candidate {key}; rerender and validate before publication"
            )
    if set(saved["parameters"]) != set(renderer.initial):
        raise ValueError("Candidate parameter surface differs from current renderer")
    seconds = saved["duration_seconds"]
    expected = renderer.render(saved["parameters"], seconds)
    fit = json.loads((directory / "candidate.fit.json").read_text(encoding="utf8"))
    check_saved_snapshot(fit, saved)
    restored = renderer.decode(
        renderer.request(command="renderSnapshot", fit=fit, seconds=seconds)["pcm"]
    )
    audio = read_wav(directory / "candidate.wav").mono()
    if audio.sample_rate != renderer.sample_rate or not np.array_equal(
        audio.samples, expected
    ):
        raise ValueError("Candidate WAV does not reproduce with the current renderer")
    if not np.array_equal(restored, expected):
        raise ValueError("Saved UI fit differs from the candidate render")
    reference = read_wav(directory / "reference.wav").mono()
    target = aligned_reference(renderer, seconds).astype(np.float32).astype(float)
    if reference.sample_rate != renderer.sample_rate or not np.array_equal(
        reference.samples, target
    ):
        raise ValueError("Report reference differs from the fixed source/alignment")
    return saved
