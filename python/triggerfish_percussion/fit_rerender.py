"""Explicit current-engine rerenders without changing the fitted reference."""

import hashlib
import json

from .audio_io import AudioBuffer, write_wav
from .workbench_fit_baseline import check_reference


def load_checked_start(renderer, path):
    """Accept engine revisions, never reference/event or parameter migrations."""
    payload = path.read_bytes()
    saved = json.loads(payload)
    previous = saved["metadata"]
    check_reference(previous, renderer.metadata)
    for key in ("recipeKey", "descriptors"):
        if previous.get(key) != renderer.metadata.get(key):
            raise ValueError(f"Fitting start {key} differs from current renderer")
    if set(saved["parameters"]) != set(renderer.initial):
        raise ValueError("Fitting start parameter surface differs from renderer")
    provenance = dict(
        path=str(path),
        sha256=hashlib.sha256(payload).hexdigest(),
        metadata=previous,
        parameters=saved["parameters"],
    )
    return saved, provenance


def prepare_current_render(renderer, source, output):
    """Validate before rendering/writing; retain complete original provenance."""
    saved, provenance = load_checked_start(renderer, source / "search.json")
    saved["rerendered_from"] = provenance
    saved["metadata"] = renderer.metadata
    audio = renderer.render(saved["parameters"], saved["duration_seconds"])
    fit = renderer.request(
        command="snapshot",
        parameters=saved["parameters"],
        name="Kick — oak medium (audition)",
    )["fit"]
    output.mkdir(parents=True, exist_ok=True)
    (output / "search.json").write_text(json.dumps(saved, indent=2), encoding="utf8")
    (output / "candidate.fit.json").write_text(
        json.dumps(fit, indent=2), encoding="utf8"
    )
    write_wav(output / "candidate.wav", AudioBuffer(audio, renderer.sample_rate))
