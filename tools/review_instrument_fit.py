"""Independent full-tail and repeated-strike audit, without publishing a fit."""

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfilt

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.reference_onset_audit import audit_onset
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

BANDS = ((40, 250), (250, 1000), (1000, 4000), (4000, 16000))
REGIONS = (
    (0, 0.03),
    (0.03, 0.12),
    (0.12, 0.5),
    (0.5, 1.5),
    (1.5, 3),
    (3, 6),
    (6, 12),
    (12, 30),
)


def band_audit(reference, candidate, rate):
    """Absolute band energy errors; floors are common and reference-only."""
    ref, synth = [], []
    for low, high in BANDS:
        sos = butter(
            2, [low, min(high, rate * 0.49)], fs=rate, btype="bandpass", output="sos"
        )
        ref.append(sosfilt(sos, reference) ** 2)
        synth.append(sosfilt(sos, candidate) ** 2)
    ref, synth = np.array(ref), np.array(synth)
    floor = max(float(ref.max()) * 1e-7, 1e-20)
    rows = []
    for start, end in REGIONS:
        first, last = round(start * rate), min(round(end * rate), len(reference))
        if last <= first:
            continue
        a = 10 * np.log10(np.maximum(ref[:, first:last].mean(axis=1), floor))
        b = 10 * np.log10(np.maximum(synth[:, first:last].mean(axis=1), floor))
        rows.append(
            dict(
                seconds=[start, last / rate],
                reference_db=a.tolist(),
                candidate_db=b.tolist(),
                difference_db=(b - a).tolist(),
            )
        )
    return dict(bands_hz=BANDS, regions=rows)


def review(target, directory):
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, directory)
        rate = renderer.sample_rate
        seconds = min(
            30.0,
            max(saved["duration_seconds"], renderer.metadata["reference"]["duration"]),
        )
        reference = aligned_reference(renderer, seconds)
        candidate = renderer.render(saved["parameters"], seconds)
        sequences = []
        for strength in (renderer.metadata["event"]["strength"], 1.0):
            hits = [
                dict(
                    time=i * 0.25,
                    strength=strength,
                    seed=(renderer.metadata["event"]["seed"] + i) & 0xFFFFFFFF,
                )
                for i in range(8)
            ]
            samples = renderer.decode(
                renderer.request(
                    command="renderSequence",
                    parameters=saved["parameters"],
                    seconds=6,
                    hits=hits,
                )["pcm"]
            )
            if not np.isfinite(samples).all():
                raise ValueError("Non-finite repeated-strike output")
            sequences.append(
                dict(
                    strength=strength,
                    peak=float(np.max(np.abs(samples))),
                    rms=float(np.sqrt(np.mean(samples**2))),
                )
            )
        report = dict(
            target=target,
            snapshot_reload_exact=True,
            renderer_sha256=renderer.metadata["rendererSha256"],
            full_duration_seconds=seconds,
            reference_onset=audit_onset(
                renderer.reference,
                rate,
                renderer.metadata["reference"]["cell"]["onset_seconds"],
            ),
            bands=band_audit(reference, candidate, rate),
            reference_peak=float(np.max(np.abs(reference))),
            candidate_peak=float(np.max(np.abs(candidate))),
            repeated_hits=sequences,
            listening_approved=False,
        )
        for label, samples in (
            ("reference-full", reference),
            ("candidate-full", candidate),
        ):
            write_wav(directory / (label + ".wav"), AudioBuffer(samples, rate))
        (directory / "full-review.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(json.dumps(report), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    review(sys.argv[1], Path(sys.argv[2]))
