"""Validate texture edits across strikes; preserve the user's modal editing.

No automatic publication. Keep absolute spectra, decay, periodicity and the
requested front protection separate rather than accepting a single scalar.
"""

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import numpy as np
import torch
from triggerfish_percussion.audio_io import AudioBuffer, read_wav, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.modulation_signature import (
    modulation_signature,
    excess_motion,
)
from fit_stretched_gong import checkpoint


def alternatives(target, base):
    yield "before", base
    if target == "crash":
        yield "clearer-packets", dict(
            base, field_phase_tilt=-1.5, field_packet_spread=2
        )
        yield "blur-only", dict(base, field_phase_tilt=-0.75)
    else:
        for depth in (0.3, 0.5):
            for blur in (base["field_phase_bandwidth"], 0.012, 0.02):
                yield f"depth-{depth}-blur-{blur:.3f}", dict(
                    base,
                    field_beat_depth=depth,
                    field_beat_rate_tilt=0.25,
                    field_phase_bandwidth=blur,
                )


def main(args):
    torch.set_num_threads(1)
    out = Path(f"build/{args.target}-texture-september")
    screen = json.loads((out / "screen.json").read_text(encoding="utf8"))
    base = screen["variants"][0]["parameters"]
    r = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        ref = aligned_reference(r, 6)
        mel = ReferenceFloorMel(ref, r.sample_rate, 60)
        texture = ModalTextureLoss(ref, r.sample_rate)
        decay = BandDecayShapeLoss(ref, r.sample_rate)
        motion = modulation_signature(ref, r.sample_rate)
        seed = r.metadata["event"]["seed"]
        rows = []
        for name, p in alternatives(args.target, base):
            seeds = []
            for offset in (0, 307, 911, 1601):
                audio = r.render(p, 6, seed + offset)
                signature = modulation_signature(audio, r.sample_rate)
                seeds.append(
                    dict(
                        seed=seed + offset,
                        mel=mel.score(audio),
                        texture=texture.score(audio),
                        decay=decay.diagnostics(audio)["shape_error_db"],
                        motion=excess_motion(signature, motion),
                        common_line_strength=signature["common_line_strength"],
                        depth=[b["depth"] for b in signature["bands"]],
                    )
                )
                if offset == 0:
                    write_wav(out / f"{name}.wav", AudioBuffer(audio, r.sample_rate))
            row = dict(
                name=name,
                parameters=p,
                seeds=seeds,
                means={
                    k: float(np.mean([s[k] for s in seeds]))
                    for k in (
                        "mel",
                        "texture",
                        "decay",
                        "motion",
                        "common_line_strength",
                    )
                },
            )
            rows.append(row)
            print(json.dumps(dict(name=name, **row["means"])), flush=True)
        (out / "validation.json").write_text(
            json.dumps(dict(reference=motion, rows=rows), indent=2), encoding="utf8"
        )
        if args.choose:
            chosen = next(row for row in rows if row["name"] == args.choose)
            loss = SpectralBloomLoss(ref, r.sample_rate)
            checkpoint(
                r, loss, out / "before", "Before texture refinement", base, ref, []
            )
            checkpoint(
                r,
                loss,
                out / "candidate",
                args.target.title() + " — texture refined",
                chosen["parameters"],
                ref,
                [
                    dict(
                        stage="texture-only refinement",
                        validation="validation.json",
                        selection=args.choose,
                    )
                ],
            )
            verify_candidate(r, out / "candidate")
            if args.target == "crash":
                original = json.loads(
                    (out / "user/user-original.fit.json").read_text(encoding="utf8")
                )
                updated = deepcopy(original)
                for node in updated["instrument"]["nodes"]:
                    for key in node["parameters"]:
                        node["parameters"][key] = chosen["parameters"][key]
                    if node["id"] == "body":
                        node["parameters"]["field_phase_tilt"] = chosen["parameters"][
                            "field_phase_tilt"
                        ]
                sounds = [
                    r.decode(
                        r.request(command="renderSnapshot", fit=f, seconds=6)["pcm"]
                    )
                    for f in (original, updated)
                ]
                front = SpectralBloomLoss(sounds[0], r.sample_rate)
                delta = (front.db(front.power(sounds[1])) - front.target)[
                    front.active, :5
                ]
                write_wav(
                    out / "candidate/user-gesture.wav",
                    AudioBuffer(sounds[1], r.sample_rate),
                )
                (out / "candidate/user-gesture.json").write_text(
                    json.dumps(
                        dict(
                            event=original["controls"]["event"],
                            front_rms_db=float(np.sqrt(np.mean(delta**2))),
                            peak_db=[
                                float(20 * np.log10(abs(s).max())) for s in sounds
                            ],
                        ),
                        indent=2,
                    )
                )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--choose")
    main(parser.parse_args())
