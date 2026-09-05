"""Keep the original audition baseline immutable across fitter/DSP revisions."""

import json

from .audio_io import AudioBuffer, read_wav, write_wav


def check_reference(previous, current):
    for key in ("sha256", "sampleRate", "referenceGainDb"):
        if previous["reference"].get(key) != current["reference"].get(key):
            raise ValueError(f"Reference {key} changed; use a new fitting directory")
    if previous["event"] != current["event"]:
        raise ValueError("Performance event changed; use a new fitting directory")
    onsets = [
        item["reference"].get("cell", {}).get("onset_seconds", 0)
        for item in (previous, current)
    ]
    if onsets[0] != onsets[1]:
        raise ValueError("Reference onset changed; use a new fitting directory")


def original_baseline(renderer, output, seconds):
    path = output / "baseline.wav"
    metadata = output / "baseline-metadata.json"
    audit = output / "audit.json"
    if path.exists():
        if metadata.exists():
            previous = json.loads(metadata.read_text())
        elif audit.exists():
            previous = json.loads(audit.read_text())["metadata"]
            metadata.write_text(json.dumps(previous, indent=2))
        else:
            raise ValueError("Existing baseline has no provenance")
        check_reference(previous, renderer.metadata)
        audio = read_wav(path).mono()
        if audio.sample_rate != renderer.sample_rate or len(audio.samples) != round(
            seconds * renderer.sample_rate
        ):
            raise ValueError(
                "Baseline duration/sample rate changed; use a new directory"
            )
        return audio.samples
    audio = renderer.render(renderer.initial, seconds)
    write_wav(path, AudioBuffer(audio, renderer.sample_rate))
    metadata.write_text(json.dumps(renderer.metadata, indent=2))
    return audio
