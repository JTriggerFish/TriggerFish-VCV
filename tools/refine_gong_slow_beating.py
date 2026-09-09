"""Test slower gong pairs without fitting pitches, levels, damping or noise.

Stage one screens rate/depth/tilt at the reference strike. Stage two checks a
named candidate across four seeds. Modulation depth by speed is kept separate
from fixed-level spectral and decay errors: no single score certifies a fit.
All rendering uses the actual workbench WASM and standard reference gesture.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.modulation_signature import modulation_signature
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def variants(base, study="slow-beating"):
    yield "before", base
    if study == "wander-hz":
        for depth in (0.3, 0.65, 1.2, 2):
            for speed in (0.2, 0.5, 1, 2):
                yield f"depth-{depth}-speed-{speed}", dict(
                    base, field_wander_hz=depth, field_wander_rate=speed
                )
        return
    for rate in (0.35, 0.65, 1, 1.5):
        for depth in (0.15, 0.3, 0.5):
            for tilt in (0, 0.25):
                yield f"rate-{rate}-depth-{depth}-tilt-{tilt}", dict(
                    base,
                    field_doublet_split=rate,
                    field_beat_depth=depth,
                    field_beat_rate_tilt=tilt,
                )


def motion_error(signature, reference):
    """Absolute relative-envelope depths; weak reference bands are excluded.

    Slow bins include smooth bloom curvature, not necessarily periodic beats.
    Return both bands of motion, not an invented ideal beating frequency.
    """
    rows = [
        (a, b)
        for a, b in zip(signature["bands"], reference["bands"])
        if b["level"] > max(r["level"] for r in reference["bands"]) * 0.03
    ]
    return {
        key: float(np.mean([abs(a[key] - b[key]) for a, b in rows]))
        for key in ("slow_depth", "fast_depth", "depth")
    }


def main(args):
    torch.set_num_threads(1)
    out = Path(f"build/{args.target}-{args.study}")
    out.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        # Archive once: repeating validation after publication must keep its baseline.
        baseline = out / "baseline.json"
        if not baseline.exists():
            baseline.write_text(json.dumps(renderer.initial, indent=2), encoding="utf8")
        base = json.loads(baseline.read_text(encoding="utf8"))
        reference = aligned_reference(renderer, 6)
        target = modulation_signature(reference, renderer.sample_rate)
        mel = ReferenceFloorMel(reference, renderer.sample_rate, 60)
        decay = BandDecayShapeLoss(reference, renderer.sample_rate)
        choices = list(variants(base, args.study))
        if args.choose:
            choices = [row for row in choices if row[0] in ("before", args.choose)]
            if len(choices) != 2:
                raise ValueError("Unknown candidate")
        offsets = (0, 307, 911, 1601) if args.choose else (0,)
        rows = []
        for name, parameters in choices:
            seeds = []
            for offset in offsets:
                audio = renderer.render(
                    parameters, 6, renderer.metadata["event"]["seed"] + offset
                )
                motion = modulation_signature(audio, renderer.sample_rate)
                seeds.append(
                    dict(
                        offset=offset,
                        mel=mel.score(audio),
                        decay=decay.diagnostics(audio)["shape_error_db"],
                        motion_error=motion_error(motion, target),
                        motion=motion,
                    )
                )
                if offset == 0:
                    write_wav(
                        out / f"{name}.wav", AudioBuffer(audio, renderer.sample_rate)
                    )
            means = {
                key: float(np.mean([s[key] for s in seeds])) for key in ("mel", "decay")
            }
            means.update(
                {
                    key: float(np.mean([s["motion_error"][key] for s in seeds]))
                    for key in ("slow_depth", "fast_depth", "depth")
                }
            )
            rows.append(
                dict(name=name, parameters=parameters, seeds=seeds, means=means)
            )
            print(json.dumps(dict(name=name, **means)), flush=True)
        filename = "validation.json" if args.choose else "screen.json"
        (out / filename).write_text(
            json.dumps(dict(reference=target, rows=rows), indent=2), encoding="utf8"
        )
        if args.choose:
            loss = SpectralBloomLoss(reference, renderer.sample_rate)
            checkpoint(
                renderer,
                loss,
                out / "before",
                args.target + " before " + args.study,
                base,
                reference,
                [],
            )
            selected = rows[-1]["parameters"]
            checkpoint(
                renderer,
                loss,
                out / "candidate",
                args.target + " — " + args.study,
                selected,
                reference,
                [dict(stage=args.study, validation="validation.json")],
            )
            verify_candidate(renderer, out / "candidate")
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=["gong", "crash"], default="gong")
    parser.add_argument(
        "--study", choices=["slow-beating", "wander-hz"], default="slow-beating"
    )
    parser.add_argument("--choose")
    main(parser.parse_args())
