"""Screen explicit gong texture trials and a measured low-core correction.

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


def clean_ring_variants(base):
    for local in (0.15, 0.35, 0.65):
        for layout in (2, 3):
            for wander in (1.5, 4):
                changed = dict(
                    base,
                    field_distribution=layout,
                    field_beat_depth=0.15,
                    field_doublet_split=2,
                    field_wander_hz=wander,
                    field_wander_rate=1,
                )
                changed.update({f"resolved_turbulence_{i}": local for i in range(4)})
                yield f"local-{local}-layout-{layout}-wander-{wander}", changed


def low_core_variants(base):
    for core in ((120, 285, 350, 540), (120, 345, 374, 540)):
        for attenuation in (0, 4, 8):
            for blur, tilt in ((0.012, 0), (0.035, -0.5)):
                changed = dict(
                    base,
                    field_beat_depth=0.15,
                    field_phase_bandwidth=blur,
                    field_phase_tilt=tilt,
                )
                for i, frequency in enumerate(core):
                    changed[f"resolved_frequency_{i}"] = frequency
                    if i:
                        changed[f"resolved_level_{i}"] -= attenuation
                yield f"core-{core[1]}-cut-{attenuation}-blur-{blur}", changed


def ring_texture_variants(base):
    for spread in (3.1252443, 4.5):
        for slope in (0, 0.2, 0.4):
            for layout in (0, 2):
                yield f"spread-{spread}-slope-{slope}-layout-{layout}", dict(
                    base,
                    field_packet_spread=spread,
                    field_turbulence_slope=slope,
                    field_distribution=layout,
                    field_beat_depth=0,
                )
    for wander in (3, 6, 12):
        for speed in (0.5, 2):
            yield f"wander-{wander}-{speed}", dict(
                base, field_wander_hz=wander, field_wander_rate=speed
            )
    for blur in (0.012, 0.035, 0.08):
        for tilt in (-1, -0.5):
            yield f"blur-{blur}-tilt-{tilt}", dict(
                base, field_phase_bandwidth=blur, field_phase_tilt=tilt
            )


def ring_balance_variants(base):
    for layout in (0, 2, 3):
        for depth in (0, 0.15, 0.3):
            if layout == 0 and depth != 0:
                continue
            yield f"layout-{layout}-depth-{depth}", dict(
                base, field_distribution=layout, field_beat_depth=depth
            )
    for wander in (0, 0.15, 1.5, 3):
        for speed in (0.2, 1.5):
            yield f"wander-{wander}-{speed}", dict(
                base, field_wander_hz=wander, field_wander_rate=speed
            )
    for spread in (1, 2, 4.5):
        yield f"spread-{spread}", dict(base, field_packet_spread=spread)


def wander_hz_variants(base):
    for depth in (0.3, 0.65, 1.2, 2):
        for speed in (0.2, 0.5, 1, 2):
            yield f"depth-{depth}-speed-{speed}", dict(
                base, field_wander_hz=depth, field_wander_rate=speed
            )


def slow_beating_variants(base):
    for rate in (0.35, 0.65, 1, 1.5):
        for depth in (0.15, 0.3, 0.5):
            for tilt in (0, 0.25):
                yield f"rate-{rate}-depth-{depth}-tilt-{tilt}", dict(
                    base,
                    field_doublet_split=rate,
                    field_beat_depth=depth,
                    field_beat_rate_tilt=tilt,
                )


def variants(base, study="slow-beating"):
    """Explicit, bounded studies; every yielded row has the complete surface."""
    studies = {
        "clean-ring": clean_ring_variants,
        "low-core": low_core_variants,
        "ring-texture": ring_texture_variants,
        "ring-balance": ring_balance_variants,
        "wander-hz": wander_hz_variants,
        "slow-beating": slow_beating_variants,
    }
    if study not in studies:
        raise ValueError("Unknown study")
    yield "before", base
    yield from studies[study](base)


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
        for key in ("slow_depth", "fast_depth", "flutter_depth", "depth")
    }


def main(args):
    torch.set_num_threads(1)
    out = args.output or Path(f"build/{args.target}-{args.study}")
    out.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        source = verify_candidate(renderer, args.source) if args.source else None
        # Archive once: repeating validation after publication must keep its baseline.
        baseline = out / "baseline.json"
        if not baseline.exists():
            initial = source["parameters"] if source else renderer.initial
            baseline.write_text(json.dumps(initial, indent=2), encoding="utf8")
        base = json.loads(baseline.read_text(encoding="utf8"))
        if set(base) != set(renderer.initial):
            raise ValueError(
                "Archived baseline uses another surface; use a fresh --output directory"
            )
        if source and source["parameters"] != base:
            raise ValueError(
                "Source differs from archived baseline; use a fresh --output"
            )
        reference = aligned_reference(renderer, 6)
        write_wav(out / "reference.wav", AudioBuffer(reference, renderer.sample_rate))
        target = modulation_signature(reference, renderer.sample_rate)
        mel = ReferenceFloorMel(reference, renderer.sample_rate, 60)
        decay = BandDecayShapeLoss(reference, renderer.sample_rate)
        shape = SpectralBloomLoss(reference, renderer.sample_rate)
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
                        envelope=float(np.linalg.norm(shape.residual(audio))),
                        motion_error=motion_error(motion, target),
                        motion=motion,
                    )
                )
                if offset == 0:
                    write_wav(
                        out / f"{name}.wav", AudioBuffer(audio, renderer.sample_rate)
                    )
            means = {
                key: float(np.mean([s[key] for s in seeds]))
                for key in ("mel", "decay", "envelope")
            }
            means.update(
                {
                    key: float(np.mean([s["motion_error"][key] for s in seeds]))
                    for key in ("slow_depth", "fast_depth", "flutter_depth", "depth")
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
            chosen = checkpoint(
                renderer,
                loss,
                out / "candidate",
                args.target + " — " + args.study,
                selected,
                reference,
                (source["history"] if source else [])
                + [
                    dict(
                        stage=args.study,
                        validation="validation.json",
                        parent=str(args.source) if args.source else "current preset",
                        selected=args.choose,
                        screen_seed=renderer.metadata["event"]["seed"],
                        validation_offsets=offsets,
                    )
                ],
            )
            if source:
                chosen.seeds = tuple(source["training_seeds"])
                chosen.save()
            verify_candidate(renderer, out / "candidate")
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=["gong", "crash"], default="gong")
    parser.add_argument(
        "--study",
        choices=[
            "slow-beating",
            "wander-hz",
            "ring-balance",
            "ring-texture",
            "low-core",
            "clean-ring",
        ],
        default="slow-beating",
    )
    parser.add_argument("--choose")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source", type=Path, help="Verified checkpoint to start from")
    main(parser.parse_args())
