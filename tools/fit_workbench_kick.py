"""Review one acoustic kick through the exact workbench signal path."""

import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from triggerfish_percussion.audio_io import AudioBuffer, read_wav, write_wav
from triggerfish_percussion.transforms import StftConfig, stft
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.short_drum_fit_loss import ShortDrumLoss
from triggerfish_percussion.workbench_fit_baseline import (
    original_baseline,
    check_reference,
)
from workbench_fit_report import write_report
from kick_fit_validation import self_test, validate_candidate
from kick_fit_stages import refine
from kick_fit_start import load_start
from kick_playability import check_playability
from kick_fit_sources import source_audit
from triggerfish_percussion.modal_fit_initialization import (
    spectral_mode_candidates,
    reference_modal_starts,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "build/workbench-wasm/site/kick-review"
SECONDS = 1.2


def inspect(samples, rate):
    transform = stft(samples, rate, StftConfig(8192, 256))
    power = transform.power
    rows = []
    for start, end in ((0, 0.03), (0.03, 0.1), (0.1, 0.25), (0.25, 0.6), (0.6, 1.2)):
        spectrum = power[
            :, (transform.times_seconds >= start) & (transform.times_seconds < end)
        ].mean(axis=1)
        peaks, _ = find_peaks(spectrum)
        peaks = sorted(peaks, key=lambda i: spectrum[i], reverse=True)[:5]
        signal = samples[round(start * rate) : round(end * rate)]
        rows.append(
            dict(
                seconds=[start, end],
                rms_db=float(10 * np.log10(max(np.mean(signal**2), 1e-30))),
                peaks_hz=transform.frequencies_hz[peaks].tolist(),
            )
        )
    return rows


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", ROOT)
    try:
        if os.environ.get("TF_KICK_SELF_TEST") == "1":
            self_test(renderer, OUTPUT / "synthetic-recovery")
            return
        if os.environ.get("TF_KICK_AUDIT_ONLY") == "1":
            write_report(OUTPUT, kind="kick", renderer=renderer)
            return
        if os.environ.get("TF_KICK_VERIFY_PRESET") == "1":
            saved = json.loads((OUTPUT / "search.json").read_text())
            check_reference(saved["metadata"], renderer.metadata)
            expected = read_wav(OUTPUT / "candidate.wav").mono().samples
            actual = renderer.render(renderer.initial, SECONDS)
            if not np.array_equal(actual, expected):
                raise RuntimeError(
                    "Reference-target preset differs from the audition render"
                )
            print(
                "Reference-target preset and audition are sample-identical", flush=True
            )
            return
        rate = renderer.sample_rate
        onset = round(
            renderer.metadata["reference"]["cell"].get("onset_seconds", 0) * rate
        )
        reference = renderer.reference[onset : onset + round(SECONDS * rate)]
        reference = np.pad(
            reference, (0, max(0, round(SECONDS * rate) - len(reference)))
        )
        baseline = original_baseline(renderer, OUTPUT, SECONDS)
        write_wav(OUTPUT / "reference.wav", AudioBuffer(reference, rate))
        report = dict(
            metadata=renderer.metadata,
            parameters=renderer.initial,
            reference=inspect(reference, rate),
            baseline=inspect(baseline, rate),
        )
        (OUTPUT / "audit.json").write_text(json.dumps(report, indent=2))
        print(
            json.dumps(
                {key: value for key, value in report.items() if key != "metadata"}
            ),
            flush=True,
        )
        loss = ShortDrumLoss(reference, rate)
        search = Search(
            renderer,
            loss,
            OUTPUT,
            seconds=SECONDS,
            name="Acoustic kick",
            seeds=(None, renderer.metadata["event"]["seed"] + 11),
        )
        resume = os.environ.get("TF_KICK_FIT_RESUME")
        if resume:
            load_start(search, resume)
        elif os.environ.get("TF_KICK_SOURCE_AUDIT") != "1":
            proposals = spectral_mode_candidates(reference, rate)
            (OUTPUT / "modal-proposals.json").write_text(
                json.dumps(proposals, indent=2)
            )
            search.screen_candidates(
                "reference-informed modal placement",
                reference_modal_starts(search.parameters, proposals),
            )
        if os.environ.get("TF_KICK_SOURCE_AUDIT") == "1":
            source_audit(renderer, search.parameters, OUTPUT / "sources")
            return
        print(json.dumps(dict(baseline=loss.diagnostics(baseline))), flush=True)
        # No output-EQ variables or hidden modal templates participate in fitting.
        fine = os.environ.get("TF_KICK_FINE_ONLY") == "1"
        refine(
            search,
            joint_only=fine or os.environ.get("TF_KICK_JOINT_ONLY") == "1",
            fine_only=fine,
        )
        search.save()
        restored = validate_candidate(renderer, search, loss, OUTPUT, SECONDS)
        check_playability(renderer, search.parameters, OUTPUT)
        source_audit(renderer, search.parameters, OUTPUT / "sources")
        write_report(OUTPUT, kind="kick", renderer=renderer)
        print(
            json.dumps(
                dict(candidate=inspect(restored, rate), round_trip="sample-identical")
            ),
            flush=True,
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
