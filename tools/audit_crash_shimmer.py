"""Validate a crash candidate on unseen phase seeds, dynamics and repeated hits."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.signal import periodogram

from crash_texture_diagnostics import CrashBalance, ridge_contrast, plots, spectrograms
from crash_refinement_contracts import audit_seeds
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.modulation_signature import modulation_signature
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def dominant_low(audio, rate):
    f, p = periodogram(audio[round(0.2 * rate) : round(2 * rate)], rate, window="hann")
    keep = (f >= 80) & (f < 180)
    return float(f[keep][np.argmax(p[keep])])


def run(directory, user_snapshot=None):
    torch.set_num_threads(1)
    fits = {
        name: json.loads((directory / filename).read_text(encoding="utf8"))
        for name, filename in [
            ("before", "source.fit.json"),
            ("candidate", "candidate.fit.json"),
        ]
    }
    seeds = audit_seeds(*(fit["controls"]["event"]["seed"] for fit in fits.values()))
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        reference = aligned_reference(renderer, 6)
        objective = CrashBalance(reference, renderer.sample_rate)
        signals = {"Reference": reference}
        write_wav(
            directory / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        results = {
            "objective": objective.specification,
            "audit_seeds": seeds,
            "renderer_sha256": renderer.metadata["rendererSha256"],
            "reference_ridges": ridge_contrast(reference, renderer.sample_rate),
            "reference_motion": modulation_signature(reference, renderer.sample_rate),
            "reference_low_peak_hz": dominant_low(reference, renderer.sample_rate),
        }
        user_event = None
        if user_snapshot is not None:
            user_event = json.loads(user_snapshot.read_text(encoding="utf8"))[
                "controls"
            ]["event"]
            results["user_gesture"] = dict(
                event=user_event,
                note="Different gesture from the reference; not a reference-fit score",
            )
        for name, fit in fits.items():
            saved = SavedFitRenderer(renderer, fit)
            if user_event is not None:
                gesture_audio = saved.render(saved.initial, 6, event=user_event)
                write_wav(
                    directory / f"{name}-user-gesture.wav",
                    AudioBuffer(gesture_audio, renderer.sample_rate),
                )
                results["user_gesture"][name] = dict(
                    peak_dbfs=float(
                        20 * np.log10(max(np.max(np.abs(gesture_audio)), 1e-15))
                    ),
                    energy=float(np.sum(gesture_audio**2)),
                )
            canonical = saved.render(saved.initial, 6)
            signals["Before" if name == "before" else "Candidate"] = canonical
            write_wav(
                directory / f"{name}.wav", AudioBuffer(canonical, renderer.sample_rate)
            )
            rows = []
            for seed in seeds:
                audio = saved.render(saved.initial, 6, seed)
                rows.append(
                    dict(
                        seed=seed,
                        components=objective.components(audio),
                        low_peak_hz=dominant_low(audio, renderer.sample_rate),
                        motion=modulation_signature(audio, renderer.sample_rate),
                        ridges=ridge_contrast(audio, renderer.sample_rate),
                    )
                )
            results[name] = rows
            if name != "candidate":
                continue
            for strength in [0.3, 0.5669291338582677, 0.9]:
                audio = saved.render(saved.initial, 6, event={"strength": strength})
                write_wav(
                    directory / f"velocity-{strength:.2f}.wav",
                    AudioBuffer(audio, renderer.sample_rate),
                )
            hits = [dict(time=0.5 * i) for i in range(5)]
            audio = renderer.decode(
                renderer.request(
                    command="renderSequence",
                    fit=fit,
                    seconds=6,
                    hits=hits,
                )["pcm"]
            )
            write_wav(
                directory / "repeated.wav", AudioBuffer(audio, renderer.sample_rate)
            )
            results["repeat_event"] = fit["controls"]["event"]
            results["repeat_peak_dbfs"] = float(
                20 * np.log10(max(1e-15, max(abs(audio))))
            )
            results["repeat_peak_after_master_minus12_dbfs"] = (
                results["repeat_peak_dbfs"] - 12
            )
        (directory / "audit.json").write_text(
            json.dumps(results, indent=2), encoding="utf8"
        )
        plots(signals, renderer.sample_rate, directory)
        spectrograms(reference, signals["Candidate"], renderer.sample_rate, directory)
        print(
            json.dumps(
                {
                    name: [r["components"] for r in results[name]]
                    for name in ["before", "candidate"]
                },
                indent=2,
            )
        )
        print(json.dumps({k: v for k, v in results.items() if k.startswith("repeat_")}))
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--user-snapshot", type=Path)
    args = parser.parse_args()
    run(args.directory, args.user_snapshot)
