"""Fail-closed report publication using the current exact renderer, not file age."""

import json
import numpy as np
from .audio_io import read_wav
from .workbench_fit_baseline import check_reference


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
    onset = round(
        renderer.metadata["reference"].get("cell", {}).get("onset_seconds", 0)
        * renderer.sample_rate
    )
    target = renderer.reference[onset : onset + len(expected)]
    target = (
        np.pad(target, (0, len(expected) - len(target)))
        .astype(np.float32)
        .astype(float)
    )
    if reference.sample_rate != renderer.sample_rate or not np.array_equal(
        reference.samples, target
    ):
        raise ValueError("Report reference differs from the fixed source/alignment")
    return saved
