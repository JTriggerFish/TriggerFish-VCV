"""Fixed-level, multi-seed and repeated-hit audit of an unpublished checkpoint."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.attack_ridge_loss import AttackRidgeLoss
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.low_mode_beating import LowModeBeating
from triggerfish_percussion.modulation_signature import modulation_signature


def review(target, directory, baseline_path, offsets, jtfs=False, floor_audit=False):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{target}-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, directory)
        source = json.loads(baseline_path.read_text(encoding="utf8"))
        if "instrument" in source:
            baseline = {
                key: value
                for node in source["instrument"]["nodes"]
                for key, value in node["parameters"].items()
            }
        else:
            baseline = source.get("baseline_parameters", source.get("parameters"))
        if not isinstance(baseline, dict) or set(baseline) != set(renderer.initial):
            raise ValueError("Baseline must contain the complete parameter surface")
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        rate = renderer.sample_rate
        mel = AuralossMel(reference, rate)
        attack_frames = round(0.4 * rate)
        attack_mel = AuralossMel(reference[:attack_frames], rate)
        bounded_mel = (
            {db: ReferenceFloorMel(reference, rate, db) for db in (50, 60, 70)}
            if floor_audit
            else {}
        )
        attack_ridges = AttackRidgeLoss(reference, rate)
        balance = MetallicBalanceLoss(reference, rate, "erb", True)
        decay = BandDecayShapeLoss(reference, rate)
        texture = ModalTextureLoss(reference, rate)
        low_beating = LowModeBeating(reference, rate) if target == "crash" else None
        scattering = None
        if jtfs:
            from triggerfish_percussion.perceptual_fit_losses import JtfsLoss

            scattering = JtfsLoss(reference, rate)
        seed = renderer.metadata["event"]["seed"]
        rows = []
        for offset in (0, *offsets):
            actual = (seed + offset) & 0xFFFFFFFF
            # An extra seed may have been used during fitting. Do not label it
            # held-out validation without checking the run's training provenance.
            row = dict(seed=actual, role="standard" if not offset else "additional")
            for label, parameters in (
                ("baseline", baseline),
                ("candidate", saved["parameters"]),
            ):
                samples = renderer.render(parameters, seconds, actual)
                row[label] = dict(
                    mel=mel.score(samples),
                    attack_mel=attack_mel.score(samples[:attack_frames]),
                    attack_ridge_mrstft=attack_ridges.score(samples),
                    balance=balance.diagnostics(samples),
                    decay=decay.diagnostics(samples),
                    texture=texture.score(samples),
                    peak_db=float(20 * np.log10(max(1e-15, abs(samples).max()))),
                )
                if seconds >= 6:
                    row[label]["modulation"] = modulation_signature(samples, rate)
                if scattering is not None:
                    row[label]["jtfs_below_8khz"] = scattering.score(samples)
                if low_beating is not None:
                    bands = low_beating.analyze(samples)
                    row[label]["low_beating"] = dict(
                        score=low_beating.score_rows(bands),
                        bands=[
                            {
                                k: band[k]
                                for k in (
                                    "band_hz",
                                    "power",
                                    "dominant_hz",
                                    "fast_fraction",
                                    "band_rms",
                                )
                            }
                            for band in bands
                        ],
                    )
                if bounded_mel:
                    row[label]["reference_floor_mel"] = {
                        db: loss.score(samples) for db, loss in bounded_mel.items()
                    }
            rows.append(row)
        repeated = {}
        for label, interval, strength in (
            ("quarters", 0.5, None),
            ("rapid-hard", 0.125, 1),
        ):
            hits = [
                dict(time=index * interval, seed=(seed + index) & 0xFFFFFFFF)
                for index in range(8)
            ]
            if strength is not None:
                hits = [dict(hit, strength=strength) for hit in hits]
            samples = renderer.decode(
                renderer.request(
                    command="renderSequence",
                    parameters=saved["parameters"],
                    seconds=12,
                    hits=hits,
                )["pcm"]
            )
            if not np.isfinite(samples).all():
                raise ValueError("Nonfinite repeated-hit output")
            write_wav(directory / f"{label}.wav", AudioBuffer(samples, rate))
            repeated[label] = dict(
                peak_db=float(20 * np.log10(max(1e-15, abs(samples).max()))),
                hit_window_rms_db=[
                    float(
                        10
                        * np.log10(
                            max(
                                1e-30,
                                np.mean(
                                    samples[
                                        round(hit["time"] * rate) : round(
                                            (hit["time"] + interval) * rate
                                        )
                                    ]
                                    ** 2
                                ),
                            )
                        )
                    )
                    for hit in hits
                ],
            )
        report = dict(
            baseline=str(baseline_path),
            baseline_parameters=baseline,
            reference=renderer.metadata["reference"],
            attack_ridge_objective=attack_ridges.specification,
            texture_objective=texture.specification,
            low_beating_diagnostic=low_beating.specification if low_beating else None,
            reference_low_beating=(
                [
                    {
                        k: band[k]
                        for k in (
                            "band_hz",
                            "power",
                            "dominant_hz",
                            "fast_fraction",
                            "band_rms",
                        )
                    }
                    for band in low_beating.target
                ]
                if low_beating
                else None
            ),
            scattering_objective=scattering.specification if scattering else None,
            reference_floor_objectives={
                db: loss.specification for db, loss in bounded_mel.items()
            },
            seeds=rows,
            repeated=repeated,
            listening_approved=False,
        )
        (directory / "refit-audit.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(dict(seeds=rows, repeated=repeated)), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("directory", type=Path)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--seed-offsets", nargs=3, type=int, required=True)
    parser.add_argument(
        "--jtfs",
        action="store_true",
        help="Optional independent WaveSpin audit (16-kHz analysis, excludes highs above 8 kHz)",
    )
    parser.add_argument(
        "--floor-audit",
        action="store_true",
        help="Compare reference-fixed Mel floors at 50, 60 and 70 dB",
    )
    args = parser.parse_args()
    if len(set(args.seed_offsets)) != 3 or any(
        not 0 < n < 2**32 for n in args.seed_offsets
    ):
        parser.error("Seed offsets must be distinct, nonzero 32-bit integers")
    review(
        args.target,
        args.directory,
        args.baseline,
        args.seed_offsets,
        args.jtfs,
        args.floor_audit,
    )
