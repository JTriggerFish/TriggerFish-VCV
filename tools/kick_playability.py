"""Offline retrigger/velocity checks through the same block API as live audio."""

import json

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav


def check_playability(renderer, parameters, output):
    rate = renderer.sample_rate
    observations = []
    for interval, strength in ((0.5, 0.5), (0.125, 0.5), (0.125, 1)):
        name = f"hits-{interval:g}s-strength-{strength:g}"
        hits = [
            dict(time=i * interval, strength=strength, seed=1449 + i) for i in range(8)
        ]
        seconds = hits[-1]["time"] + 1.2
        audio = renderer.decode(
            renderer.request(
                command="renderSequence",
                parameters=parameters,
                seconds=seconds,
                hits=hits,
            )["pcm"]
        )
        if not np.all(np.isfinite(audio)):
            raise RuntimeError("Non-finite retrigger output")
        write_wav(output / f"{name}.wav", AudioBuffer(audio, rate))
        attacks = [
            audio[round(hit["time"] * rate) : round((hit["time"] + 0.03) * rate)]
            for hit in hits
        ]
        observations.append(
            dict(
                name=name,
                peak=float(np.max(np.abs(audio))),
                attack_rms_db=[
                    float(10 * np.log10(max(np.mean(a * a), 1e-30))) for a in attacks
                ],
            )
        )
    # The exact same first event rendered in a sequence must equal a single hit;
    # this catches a hidden reset, different gain path, or incorrect block offset.
    sequence = renderer.decode(
        renderer.request(
            command="renderSequence",
            parameters=parameters,
            seconds=1.2,
            hits=[dict(time=0)],
        )["pcm"]
    )
    if not np.array_equal(sequence, renderer.render(parameters, 1.2)):
        raise RuntimeError("Sequence and single-hit render paths disagree")
    (output / "playability.json").write_text(json.dumps(observations, indent=2))
    check_resonance_mix(renderer, parameters, output)
    return observations


def check_resonance_mix(renderer, parameters, output):
    """Prominence must scale only the resonator observation, without auto-level."""
    dry = renderer.render(dict(parameters, resonance_level=0), 1.2)
    wet = renderer.render(dict(parameters, resonance_level=1), 1.2) - dry
    results = []
    for level in (0, 0.25, 0.5, 1, 2):
        audio = renderer.render(dict(parameters, resonance_level=level), 1.2)
        error = float(np.max(np.abs(audio - (dry + level * wet))))
        if error > 2e-4 * max(1.0, float(np.max(np.abs(audio)))):
            raise RuntimeError(
                "Resonance prominence changed more than its observation gain"
            )
        results.append(dict(resonance_level=level, affine_error=error))
        write_wav(
            output / f"resonance-{level:g}.wav",
            AudioBuffer(audio, renderer.sample_rate),
        )
    (output / "resonance-mix.json").write_text(json.dumps(results, indent=2))
