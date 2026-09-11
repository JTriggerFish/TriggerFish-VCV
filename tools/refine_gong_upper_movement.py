"""Coordinate search for upper metallic movement, guarded by layered timing.

Uses the actual WASM engine; no added noise, hidden gain, or modal placement.
Modulation statistics are ranking diagnostics, not perceptual-equivalence proof.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from audit_gong_sizzle import brightness, local_texture, ridge_contrast
from fit_stretched_gong import checkpoint

CHOICES = {
    "field_beat_rate_tilt": [0.5, 0.75, 1],
    "field_beat_depth": [0.3, 0.5],
    "field_satellite_density": [0.2, 0.7, 1],
    "field_packet_spread": [1, 1.8, 4],
    "field_motion_depth": [0.5, 1.5],
    "field_motion_rate": [40, 200],
    "field_motion_sharing": [0, 0.75],
    "field_phase_bandwidth": [0.01, 0.075],
}


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = json.loads(args.source.read_text())["parameters"]
        ref = aligned_reference(r, 6)
        loss = LayeredBandLoss(ref, r.sample_rate)
        texture = ModalTextureLoss(ref, r.sample_rate)
        target_mod = np.maximum(local_texture(ref, texture)[1], 1e-12)
        target_centroid = brightness(ref, r.sample_rate)[0]["centroid_hz"]
        target_ridge = np.array(ridge_contrast(ref, r.sample_rate))
        args.output.mkdir(parents=True, exist_ok=True)
        rows = []

        def measure(p):
            values = []
            for seed in (1675, 1982):
                audio = r.render(p, 6, seed)
                db = loss.envelopes(audio)
                mod = np.maximum(local_texture(audio, texture)[1], 1e-12)
                centroid = brightness(audio, r.sample_rate)[0]["centroid_hz"]
                values.append(
                    dict(
                        seed=seed,
                        timing=loss.score_db(db),
                        shape=loss.score_db(db, True),
                        low_db=db[:3].tolist(),
                        ridge=ridge_contrast(audio, r.sample_rate),
                        modulation_db=float(
                            np.sqrt(np.mean((10 * np.log10(mod / target_mod)) ** 2))
                        ),
                        modulation_power=mod.tolist(),
                        centroid_hz=centroid,
                        colour_octaves=float(abs(np.log2(centroid / target_centroid))),
                    )
                )
            return values

        original = measure(base)

        def record(name, p, values):
            eligible = True
            for v, b in zip(values, original):
                eligible &= v["timing"] <= b["timing"] + 0.5
                eligible &= v["shape"] <= b["shape"] + 0.4
                eligible &= (
                    np.sqrt(np.mean((np.array(v["low_db"]) - b["low_db"]) ** 2)) <= 1.5
                )
                eligible &= np.all(np.array(v["ridge"]) >= target_ridge - 0.9)
            scores = [
                v["timing"] + 0.5 * v["modulation_db"] + 6 * v["colour_octaves"]
                for v in values
            ]
            row = dict(
                name=name,
                parameters=p,
                metrics=values,
                eligible=bool(eligible),
                score=float(np.mean(scores) + 0.25 * max(scores)),
            )
            rows.append(row)
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            print(
                json.dumps({k: row[k] for k in ("name", "score", "eligible")}),
                flush=True,
            )
            return row

        best = record("baseline", base, original)
        for round_number in range(args.rounds):
            start = dict(best["parameters"])
            for key, values in CHOICES.items():
                for value in values:
                    if value == start[key]:
                        continue
                    p = dict(start, **{key: value})
                    row = record(f"round {round_number}: {key}={value}", p, measure(p))
                    if row["eligible"] and row["score"] < best["score"]:
                        best = row
        checkpoint(
            r,
            loss,
            args.output / "candidate",
            "Gong — shaped bloom and sizzle",
            best["parameters"],
            ref,
            [
                dict(
                    stage="guarded upper movement",
                    selected=best["name"],
                    score=best["score"],
                    baseline_score=rows[0]["score"],
                    timing=loss.specification,
                    choices=CHOICES,
                    ranking="mean plus quarter worst of timing + 0.5 modulation dB + 6 centroid octaves",
                    guards="timing +0.5 dB, shape +0.4 dB, low change RMS 1.5 dB, ridge >= reference -0.9 dB",
                )
            ],
        )
        verify_candidate(r, args.output / "candidate")
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=2)
    run(parser.parse_args())
