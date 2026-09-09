"""Reference-led texture screen; no free modal placement, decay or gain fitting."""

import argparse
import json
import os
from pathlib import Path
import numpy as np
import torch
from scipy.signal import find_peaks, periodogram, resample_poly
from math import gcd
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.modulation_signature import (
    modulation_signature,
    excess_motion,
)
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss


def variants(target, base):
    yield "before", dict(base)
    if target == "crash":
        for tilt in (-0.75, -1.5, -2):
            yield f"blur-tilt-{tilt}", dict(base, field_phase_tilt=tilt)
        for blur in (0.01, 0.02):
            yield f"less-blur-{blur}", dict(base, field_phase_bandwidth=blur)
        for spread in (1.6, 2.0):
            yield f"clearer-{spread}", dict(
                base, field_phase_tilt=-1.5, field_packet_spread=spread
            )
        for scale in (0.65, 0.8):
            yield f"clear-rings-{scale}", dict(
                base,
                field_phase_tilt=-1.5,
                field_turbulence=base["field_turbulence"] * scale,
                field_packet_spread=base["field_packet_spread"] / scale,
            )
        for level in (-42, -36):
            yield f"low-companion-{level}", dict(
                base,
                field_phase_tilt=-1.5,
                field_packet_spread=2,
                resolved_frequency_24=81,
                resolved_level_24=level,
                resolved_turbulence_24=0.1,
                resolved_allocation_24=0,
            )
    else:
        for depth in (0.15, 0.3, 0.5):
            for tilt in (0.25, 0.5):
                yield f"gentle-doublets-{depth}-{tilt}", dict(
                    base, field_beat_depth=depth, field_beat_rate_tilt=tilt
                )
        for depth in (0.03, 0.1, 0.3):
            for speed in (0.8, 2):
                yield f"wander-hz-{depth}-{speed}", dict(
                    base, field_wander_hz=depth, field_wander_rate=speed
                )
        for layout in (0, 3):
            for drift in (0, 0.1):
                yield f"layout-{layout}-drift-{drift}", dict(
                    base,
                    field_distribution=layout,
                    field_doublet_split=0.8,
                    field_beat_depth=0.15,
                    field_wander_hz=drift,
                    field_wander_rate=0.8,
                )


def low_peaks(audio, rate):
    divisor = gcd(int(rate), 2000)
    x = resample_poly(audio, 2000 // divisor, int(rate) // divisor)
    rows = []
    for start in (0.15, 2.15, 4):
        f, p = periodogram(
            x[round(start * 2000) : round((start + 2) * 2000)],
            fs=2000,
            window="hann",
            nfft=16000,
        )
        indices, _ = find_peaks(p)
        indices = [i for i in indices if 40 < f[i] < 180]
        top = sorted(indices, key=lambda i: p[i], reverse=True)[:8]
        rows.append(
            dict(
                start=start,
                window_seconds=2,
                bin_hz=0.125,
                resolving_width_hz=0.5,
                peaks=[[float(f[i]), float(10 * np.log10(p[i]))] for i in top],
            )
        )
    return rows


def main(args):
    torch.set_num_threads(1)
    out = Path(f"build/{args.target}-texture-september")
    out.mkdir(parents=True, exist_ok=True)
    r = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        base = dict(r.initial)
        if args.fit:
            saved = json.loads(args.fit.read_text(encoding="utf8"))
            base.update(
                {
                    k: v
                    for n in saved["instrument"]["nodes"]
                    for k, v in n["parameters"].items()
                }
            )
        reference = aligned_reference(r, 6)
        mel = ReferenceFloorMel(reference, r.sample_rate, 60)
        texture = ModalTextureLoss(reference, r.sample_rate)
        modulation = modulation_signature(reference, r.sample_rate)
        write_wav(out / "reference.wav", AudioBuffer(reference, r.sample_rate))
        records = []
        original = r.render(base, 6)
        front = SpectralBloomLoss(original, r.sample_rate)
        for name, p in variants(args.target, base):
            audio = r.render(p, 6)
            signature = modulation_signature(audio, r.sample_rate)
            delta = (front.db(front.power(audio)) - front.target)[front.active, :5]
            row = dict(
                name=name,
                parameters=p,
                mel=mel.score(audio),
                texture=texture.score(audio),
                motion=excess_motion(signature, modulation),
                signature=signature,
                front_db=float(np.sqrt(np.mean(delta**2))),
            )
            records.append(row)
            write_wav(out / f"{name}.wav", AudioBuffer(audio, r.sample_rate))
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in ("name", "mel", "texture", "motion", "front_db")
                    }
                ),
                flush=True,
            )
        report = dict(
            reference_signature=modulation,
            reference_low_peaks=low_peaks(reference, r.sample_rate),
            initial_low_peaks=low_peaks(original, r.sample_rate),
            metadata=r.metadata,
            variants=records,
        )
        (out / "screen.json").write_text(json.dumps(report, indent=2), encoding="utf8")
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--fit", type=Path)
    main(parser.parse_args())
