"""Screen low-packet clarity without sacrificing the established gong bloom.

Exact WASM renders, fixed reference level and gesture. The texture study keeps
all centres; low-core studies explicitly replace only the first four/five
centres with reference-guided hypotheses. No per-mode damping, extra knots,
gain matching or engine changes. Ranking uses
the existing Auraloss attack MR-STFT below 1.5 kHz; whole-sound Mel and band
envelope errors are separate acceptance guards, not folded into that ranking.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.signal import butter, find_peaks, sosfilt, welch

from triggerfish_percussion.attack_ridge_loss import AttackRidgeLoss
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def low_peaks(audio, rate):
    segment = audio[round(0.04 * rate) : round(0.5 * rate)]
    f, p = welch(segment, rate, nperseg=8192, nfft=32768)
    peaks, _ = find_peaks(p)
    peaks = [i for i in peaks if 80 < f[i] < 1000]
    strongest = sorted(peaks, key=lambda i: p[i], reverse=True)[:10]
    return [
        dict(hz=float(f[i]), db=float(10 * np.log10(max(p[i], 1e-20))))
        for i in strongest
    ]


def candidates(base, study):
    yield "before", dict(base)
    if study == "front":
        for concentration in (0.05, 0.1, 0.2, 0.4, 0.7, 1):
            for rate in (1, 3, 8, 16):
                p = dict(base, bloom_energy_acceleration=concentration, bloom_rate=rate)
                for i, f in enumerate((121, 288, 344, 374, 537)):
                    p[f"resolved_frequency_{i}"] = f
                    p[f"resolved_turbulence_{i}"] = 0.65
                for i, delta in enumerate((0, -7, 2, 2, -4)):
                    p[f"resolved_level_{i}"] = base[f"resolved_level_{i}"] + delta
                yield f"concentration-{concentration}-rate-{rate}", p
        return
    if study in ("low-core", "five-core"):
        for local in (0.25, 0.4, 0.65):
            for balance in (-2, 0, 2):
                for rate in (0.7, 1, 1.3):
                    p = dict(base, bloom_rate=base["bloom_rate"] * rate)
                    for i, frequency in enumerate((121, 344, 374, 537)):
                        p[f"resolved_frequency_{i}"] = frequency
                        p[f"resolved_turbulence_{i}"] = local
                    p["resolved_level_1"] += 4 + balance
                    p["resolved_level_2"] += -2 + balance
                    p["resolved_level_3"] -= 3
                    if study == "five-core":
                        for i, frequency in enumerate((121, 288, 344, 374, 537)):
                            p[f"resolved_frequency_{i}"] = frequency
                            p[f"resolved_turbulence_{i}"] = local
                        for i, delta in enumerate((0, -7, balance, balance, -4)):
                            p[f"resolved_level_{i}"] = (
                                base[f"resolved_level_{i}"] + delta
                            )
                    yield f"pair-local-{local}-balance-{balance}-rate-{rate}", p
        return
    for local in (0, 0.15, 0.35, 0.65):
        for lift in (0, 4, 8):
            for rate in (0.5, 1, 1.5):
                p = dict(base, bloom_rate=base["bloom_rate"] * rate)
                p["resolved_level_0"] += lift
                p.update({f"resolved_turbulence_{i}": local for i in range(4)})
                yield f"local-{local}-core-{lift}-rate-{rate}", p


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        baseline_path = args.output / "baseline.json"
        if baseline_path.exists():
            base = json.loads(baseline_path.read_text(encoding="utf8"))
        else:
            base = dict(r.initial)
            baseline_path.write_text(json.dumps(base, indent=2), encoding="utf8")
        ref = aligned_reference(r, 6)
        lowpass = butter(6, 1500, fs=r.sample_rate, output="sos")
        attack = AttackRidgeLoss(sosfilt(lowpass, ref), r.sample_rate, 0.5)
        mel = ReferenceFloorMel(ref, r.sample_rate, 60)
        bloom = SpectralBloomLoss(ref, r.sample_rate)
        rows = []
        checkpoint(
            r,
            bloom,
            args.output / "baseline",
            "Gong before attack refinement",
            base,
            ref,
            [],
        )
        print(
            json.dumps(dict(reference_peaks=low_peaks(ref, r.sample_rate))), flush=True
        )
        for name, parameters in candidates(base, args.study):
            if args.baseline_only and name != "before":
                break
            audio = r.render(parameters, 6)
            diagnostics = bloom.diagnostics(audio)
            row = dict(
                name=name,
                parameters=parameters,
                attack=attack.score(sosfilt(lowpass, audio)),
                mel=mel.score(audio),
                envelope=diagnostics["envelope_rms_db"],
                rise=diagnostics["rise_rms_db"],
                peaks=low_peaks(audio, r.sample_rate),
            )
            rows.append(row)
            print(
                json.dumps(
                    {k: v for k, v in row.items() if k not in ("parameters", "peaks")}
                ),
                flush=True,
            )
            (args.output / "screen.json").write_text(
                json.dumps(rows, indent=2), encoding="utf8"
            )
        before = rows[0]
        eligible = [
            row
            for row in rows
            if row["mel"] <= before["mel"] * 1.08
            and row["envelope"] <= before["envelope"] * 1.10
            and row["rise"] <= before["rise"] * 1.10
        ]
        for index, row in enumerate(
            sorted(eligible, key=lambda row: row["attack"])[:3]
        ):
            checkpoint(
                r,
                bloom,
                args.output / f"candidate-{index}",
                "Gong pitched attack trial",
                row["parameters"],
                ref,
                [
                    dict(
                        stage="pitched attack screen",
                        selected=row["name"],
                        attack_specification=attack.specification,
                        lowpass_hz=1500,
                        constraints=dict(
                            mel_factor=1.08, envelope_factor=1.10, rise_factor=1.10
                        ),
                    )
                ],
            )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("build/gong-pitched-attack-e2a7eb6")
    )
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument(
        "--study",
        choices=["texture", "low-core", "five-core", "front"],
        default="texture",
    )
    run(parser.parse_args())
