"""Audit the exact published workbench patch without replacing it."""

import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.band_region_audit import BandRegionAudit
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from kick_control_span import (
    control_probes,
    contact_probes,
    coverage_layout,
    observation_probes,
)
from kick_fit_sources import ROUTES

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "build/kick-diagnosis"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", ROOT)
    try:
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset : onset + round(1.2 * rate)]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        audit = BandRegionAudit(reference, rate)
        initial = dict(renderer.initial)
        rows, parts = [], []

        def measure(name, values, save=False):
            audio = renderer.render(values, 1.2)
            row = dict(name=name, parameters=values, **audit.measure(audio))
            rows.append(row)
            if save:
                write_wav(OUTPUT / f"{name}.wav", AudioBuffer(audio, rate))
            return audio

        write_wav(OUTPUT / "reference.wav", AudioBuffer(reference, rate))
        baseline = measure("published", initial, True)
        for route in ROUTES:
            isolated = dict(initial, **{key: 0 for key in ROUTES if key != route})
            parts.append(measure(route, isolated, True))
        no_noise = measure(
            "without-source-noise", dict(initial, contact_noise_level=0), True
        )
        write_wav(
            OUTPUT / "source-noise-contribution.wav",
            AudioBuffer(baseline - no_noise, rate),
        )
        for index, (name, values) in enumerate(control_probes(initial)):
            measure(name, values)
        for name, values in contact_probes(initial):
            measure(name, values)
        for name, values in observation_probes(initial):
            measure(name, values)
        covered = coverage_layout(initial)
        measure("coverage-layout", covered, True)
        for name, values in control_probes(covered):
            measure("coverage-" + name, values)
        report = dict(
            metadata=renderer.metadata,
            rows=rows,
            source_sum_max_error=float(
                np.max(np.abs(baseline - np.sum(parts, axis=0)))
            ),
        )
        (OUTPUT / "audit.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(
            json.dumps(
                dict(
                    renders=len(rows),
                    source_sum_error=report["source_sum_max_error"],
                    baseline=rows[0],
                    best=sorted(rows[5:], key=lambda r: r["worst_audible_error_db"])[
                        :3
                    ],
                )
            ),
            flush=True,
        )
        if os.environ.get("TF_KICK_SPAN_REFIT") == "1":
            from kick_span_refinement import refine_span

            refine_span(renderer, audit, reference, OUTPUT)
        if os.environ.get("TF_KICK_MEASURED_LAYOUT") == "1":
            from kick_measured_span import refine_measured

            refine_measured(renderer, audit, reference, OUTPUT)
        from kick_diagnostic_plots import plot_diagnosis

        plot_diagnosis(OUTPUT)
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
