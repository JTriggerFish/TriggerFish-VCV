"""Controlled texture alternatives around a saved series, not ridge fitting.

Use the actual WASM renderer and existing auditory-band modulation diagnostics.
All modal frequencies/prominences and output gains stay fixed. Allocation
experiments use one smooth frequency law, never individually fitted weights.
"""

import argparse
import json
import os
from pathlib import Path
from itertools import product

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from refine_edited_gong import Comparison, plot_signals


def variations(p, joint=False, finish=False):
    yield "Your gong test", dict(p)
    if finish:
        for rate, decay in product((2.8, 3.2, 3.6), (1.75, 2, 2.25)):
            yield f"Bloom {rate}, T60 {decay}", dict(
                p,
                field_motion_depth=2.2,
                field_packet_spread=1.4,
                bloom_rate=rate,
                body_decay_seconds_7=decay,
            )
        return
    if joint:
        for movement, width in product((1.5, 2.2, 2.8), (1.1, 1.4, 1.8)):
            yield f"Movement {movement}, spread {width}", dict(
                p, field_motion_depth=movement, field_packet_spread=width
            )
        return
    grid = {
        "field_motion_depth": (0, 0.75, 2.2, 2.8),
        "field_motion_rate": (80, 120),
        "field_beat_depth": (0, 0.35),
        "field_packet_spread": (1, 1.4, 2.3),
        "field_satellite_density": (0.5, 0.75),
        "field_distribution": (0, 1, 2),
    }
    for key, values in grid.items():
        for value in values:
            yield f"{key}={value}", dict(p, **{key: value})
    for slope in (0.5, 1):
        trial = dict(p)
        for i in range(32):
            trial[f"resolved_allocation_{i}"] = float(
                np.clip((p[f"resolved_frequency_{i}"] / 1000) ** slope, 0.25, 4)
            )
        yield f"Upper allocation slope {slope}", trial


def texture_features(texture, audio):
    features, _ = texture.features(audio)
    # Keep the slow bloom separate: compare 8–32 and 32–128 Hz modulation
    # during bloom/tail, across six auditory bands. Log-power -> dB.
    return 10 * features.reshape(6, 3, 4)[:, 1:, 1:3]


def run(args):
    args.output.mkdir(parents=True, exist_ok=True)
    fit = json.loads(args.fit.read_text())
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(r, fit)
        seeds = (fit["controls"]["event"]["seed"], 1982)
        reference = aligned_reference(r, 6)
        baseline = [saved.render(saved.initial, 6, seed) for seed in seeds]
        targets = [Comparison(reference, x, r.sample_rate) for x in baseline]
        texture = ModalTextureLoss(
            reference, r.sample_rate, centres=[2500, 3500, 4500, 6500, 9500, 12500]
        )
        target_features = texture_features(texture, reference)
        rows = []
        for name, p in variations(saved.initial, args.joint, args.finish):
            metrics = []
            for seed, target in zip(seeds, targets):
                audio = saved.render(p, 6, seed)
                m = target.metrics(audio)
                features = texture_features(texture, audio)
                m["texture_db"] = float(
                    np.sqrt(np.mean((features - target_features) ** 2))
                )
                m["modulation_error_db"] = (features - target_features).tolist()
                metrics.append(m)
            row = dict(
                name=name,
                parameters=p,
                metrics=metrics,
                timing=float(np.mean([m["score"] for m in metrics])),
                texture=float(np.mean([m["texture_db"] for m in metrics])),
                body=float(np.mean([m["body_envelope_db"] for m in metrics])),
            )
            rows.append(row)
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            print(
                json.dumps(
                    {k: v for k, v in row.items() if k not in ("parameters", "metrics")}
                ),
                flush=True,
            )
        write_wav(args.output / "reference.wav", AudioBuffer(reference, r.sample_rate))
        write_wav(args.output / "edited.wav", AudioBuffer(baseline[0], r.sample_rate))
        (args.output / "source.fit.json").write_text(json.dumps(fit, indent=2))
        (args.output / "provenance.json").write_text(
            json.dumps(
                dict(
                    source_id=fit["id"],
                    event=fit["controls"]["event"],
                    seeds=seeds,
                    renderer_sha256=r.metadata["rendererSha256"],
                    texture=texture.specification,
                    target=targets[0].spectral.specification,
                    interpretation="Controlled diagnostic alternatives; no scalar score proves naturalness",
                ),
                indent=2,
            )
        )
        if args.finish:
            winner = min(rows, key=lambda row: np.hypot(row["timing"], row["texture"]))
            p = winner["parameters"]
            if any(
                p[k] != v for k, v in saved.initial.items() if k.startswith("resolved_")
            ):
                raise ValueError("Refinement changed a painted mode")
            candidate = saved.snapshot(p, "Gong — upper texture refinement")
            audio = saved.render(p, 6)
            replay = r.decode(
                r.request(command="renderSnapshot", fit=candidate, seconds=6)["pcm"]
            )
            if not np.array_equal(audio, replay):
                raise ValueError("Saved fit does not reproduce the rendered candidate")
            write_wav(args.output / "candidate.wav", AudioBuffer(audio, r.sample_rate))
            (args.output / "candidate.fit.json").write_text(
                json.dumps(candidate, indent=2)
            )
            (args.output / "selected.json").write_text(json.dumps(winner, indent=2))
            plot_signals(
                {
                    "Reference": reference,
                    "Your gong test": baseline[0],
                    "Candidate": audio,
                },
                r.sample_rate,
                args.output,
            )
            print(
                json.dumps(
                    dict(
                        selected=winner["name"],
                        changed={
                            k: [saved.initial[k], v]
                            for k, v in p.items()
                            if v != saved.initial[k]
                        },
                    )
                ),
                flush=True,
            )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--joint", action="store_true")
    parser.add_argument("--finish", action="store_true")
    run(parser.parse_args())
