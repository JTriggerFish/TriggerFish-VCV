"""Read-only counterfactual audit of the loss that suppressed the gong's low body.

Uses saved actual-WASM renders and fresh one-factor renders. Never publishes a
preset. Frequency slices diagnose this case, not fitted frequency exceptions.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, welch

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.regional_spectrum_audit import RegionalSpectrumAudit
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def load(directory):
    parameters = json.loads((directory / "search.json").read_text())["parameters"]
    wav = read_wav(directory / "candidate.wav").mono()
    return parameters, wav.samples, wav.sample_rate


def pitch_slices(audio, rate):
    """Report measured early peaks, not an inferred perceptual fundamental."""
    rows = []
    for start, end in ((0.04, 0.2), (0.2, 0.5), (0.5, 1), (2, 4)):
        segment = audio[round(start * rate) : round(end * rate)]
        f, p = welch(segment, rate, nperseg=min(8192, len(segment)), nfft=32768)
        peaks, _ = find_peaks(p)
        peaks = peaks[(f[peaks] >= 60) & (f[peaks] <= 900)]
        peaks = sorted(peaks, key=lambda i: p[i], reverse=True)[:8]
        rows.append(
            dict(
                seconds=[start, end],
                peaks_hz=[float(f[i]) for i in peaks],
                peak_density_db=[float(10 * np.log10(p[i])) for i in peaks],
                band_power_db={
                    f"{lo}-{hi}": float(
                        10
                        * np.log10(
                            max(1e-20, p[(f >= lo) & (f < hi)].sum() * (f[1] - f[0]))
                        )
                    )
                    for lo, hi in ((80, 180), (180, 300), (300, 500), (500, 900))
                },
            )
        )
    return rows


def describe(audio, rate, losses, spectrum):
    return dict(
        objectives={name: loss.attribution(audio) for name, loss in losses.items()},
        pitch_slices=pitch_slices(audio, rate),
        spectrum=spectrum.measure(audio),
    )


def run(args):
    old, incoming, rate = load(args.incoming)
    fitted, rejected, fitted_rate = load(args.rejected)
    reference_wav = read_wav(args.rejected / "reference.wav").mono()
    reference = reference_wav.samples
    if rate != fitted_rate or rate != reference_wav.sample_rate:
        raise ValueError("Mismatched sample rates")
    if not np.array_equal(
        reference, read_wav(args.incoming / "reference.wav").mono().samples
    ):
        raise ValueError("Saved references differ")
    losses = dict(
        uniform=LayeredBandLoss(reference, rate),
        weighted=LayeredBandLoss(reference, rate, audibility=True),
    )
    spectrum = RegionalSpectrumAudit(
        reference, rate, ((0.04, 0.2), (0.2, 0.5), (0.5, 1), (1, 2), (2, 4), (4, 6))
    )
    rows = {
        "reference": describe(reference, rate, losses, spectrum),
        "incoming": describe(incoming, rate, losses, spectrum),
        "rejected": describe(rejected, rate, losses, spectrum),
    }
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        verify_candidate(r, args.incoming)
        verify_candidate(r, args.rejected)
        trials = {
            "only_cut_first_bar": dict(
                old, resolved_level_0=fitted["resolved_level_0"]
            ),
            "restore_first_bar_in_rejected": dict(
                fitted, resolved_level_0=old["resolved_level_0"]
            ),
            "only_replace_all_bars": dict(
                old,
                **{k: v for k, v in fitted.items() if k.startswith("resolved_level_")},
            ),
        }
        for decay in (3, 5, 7):
            trials[f"incoming_low_t60_{decay}"] = dict(old, body_decay_seconds_0=decay)
        for name, p in trials.items():
            audio = r.render(p, 6, 1675)
            rows[name] = describe(audio, rate, losses, spectrum)
            rows[name]["parameters"] = p
            print(
                json.dumps(
                    dict(
                        trial=name,
                        scores={
                            k: v["score"] for k, v in rows[name]["objectives"].items()
                        },
                    )
                ),
                flush=True,
            )
    finally:
        r.close()
    args.output.mkdir(parents=True, exist_ok=True)
    payload = dict(
        reference_sha256=hashlib.sha256(reference.tobytes()).hexdigest(),
        renderer_sha256=r.metadata["rendererSha256"],
        event=r.metadata["event"],
        specification=losses["weighted"].specification,
        sources=dict(incoming=str(args.incoming), rejected=str(args.rejected)),
        rows=rows,
    )
    (args.output / "objective-audit.json").write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--incoming", type=Path, required=True)
    parser.add_argument("--rejected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
