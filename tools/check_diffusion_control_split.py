"""Exact-WASM migration and velocity sweeps; no fitting or audio normalization."""

import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, read_wav, write_wav
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def main():
    directory = Path("build/diffusion-control-split")
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        old = json.loads((directory / "baseline.json").read_text())
        expected = dict(
            old, bloom_energy_sensitivity=2 * old["bloom_energy_acceleration"]
        )
        if renderer.initial != expected:
            raise ValueError("Migration changed more than the explicit split parameter")
        previous = read_wav(directory / "gong-before.wav").mono()
        current = renderer.render(renderer.initial, 6)
        if (
            previous.sample_rate != renderer.sample_rate
            or current.shape != previous.samples.shape
        ):
            raise ValueError("Mismatched before/after render")
        relative = np.linalg.norm(current - previous.samples) / np.linalg.norm(
            previous.samples
        )
        # This is a six-second nonlinear trajectory, not a bit-identical promise.
        if relative > 0.001:
            raise ValueError(
                f"Linked preset migration changed the waveform: {relative}"
            )
        write_wav(
            directory / "gong-after.wav", AudioBuffer(current, renderer.sample_rate)
        )
        rows = []
        for sensitivity in (0, 1, 2):
            parameters = dict(renderer.initial, bloom_energy_sensitivity=sensitivity)
            for strength in (0.25, 0.5, 1):
                audio = renderer.decode(
                    renderer.request(
                        command="renderSequence",
                        parameters=parameters,
                        seconds=6,
                        hits=[dict(time=0, strength=strength)],
                    )["pcm"]
                )
                if not np.isfinite(audio).all() or not np.any(audio):
                    raise ValueError("Velocity sweep produced invalid or silent audio")
                name = f"energy-{sensitivity}-strength-{strength}.wav"
                write_wav(directory / name, AudioBuffer(audio, renderer.sample_rate))
                rows.append(
                    dict(
                        sensitivity=sensitivity,
                        strength=strength,
                        peak_db=float(20 * np.log10(np.max(np.abs(audio)))),
                        wav=name,
                    )
                )
        report = dict(
            migration_relative_waveform_error=float(relative),
            sweeps=rows,
            metadata=renderer.metadata,
            normalization=False,
        )
        (directory / "verification.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(
            json.dumps(
                dict(migration_relative_waveform_error=float(relative), sweeps=rows)
            )
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
