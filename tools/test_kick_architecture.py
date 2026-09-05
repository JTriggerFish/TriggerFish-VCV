"""Re-fit component ablations of the current single kick, separately from presets."""

import json
import os
from pathlib import Path

import numpy as np
from triggerfish_percussion.short_drum_fit_loss import ShortDrumLoss
from triggerfish_percussion.workbench_fit_baseline import check_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from kick_fit_stages import stages_for
from kick_fit_validation import validate_candidate
from kick_study_report import finish_study

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "build/workbench-wasm/site/kick-study"
START = ROOT / "build/workbench-wasm/site/kick-review/search.json"
SECONDS = 1.2


def run_case(renderer, loss, initial, name, fixed, omitted):
    output = OUTPUT / name
    output.mkdir(parents=True, exist_ok=True)
    search = Search(
        renderer,
        loss,
        output,
        seconds=SECONDS,
        name=name,
        seeds=(None, renderer.metadata["event"]["seed"] + 11),
    )
    search.parameters = dict(initial, **fixed)
    before = loss.diagnostics(search.audio(search.parameters))
    bounds = {
        key: value
        for _, stage in stages_for(initial)
        for key, value in stage.items()
        if key not in {*fixed, *omitted}
    }
    search.stage(name, bounds, 35, difference_step=0.001)
    search.save()
    audio = validate_candidate(renderer, search, loss, output, SECONDS)
    result = dict(
        name=name,
        fixed=fixed,
        fitted=list(bounds),
        deleted_source_diagnostics=before,
        fitted_diagnostics=loss.diagnostics(audio),
        evaluations=search.evaluations,
        parameters=search.parameters,
    )
    (output / "ablation.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    saved = json.loads(START.read_text())
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", ROOT)
    try:
        check_reference(saved["metadata"], renderer.metadata)
        if set(saved["parameters"]) != set(renderer.initial):
            raise ValueError("Architecture audit needs a fit of the current kick")
        onset = round(
            renderer.metadata["reference"]["cell"].get("onset_seconds", 0)
            * renderer.sample_rate
        )
        reference = renderer.reference[
            onset : onset + round(SECONDS * renderer.sample_rate)
        ]
        reference = np.pad(
            reference,
            (0, max(0, round(SECONDS * renderer.sample_rate) - len(reference))),
        )
        if os.environ.get("TF_KICK_STUDY_AUDIT") == "1":
            finish_study(renderer, reference, OUTPUT)
            return
        loss = ShortDrumLoss(reference, renderer.sample_rate)
        initial = saved["parameters"]
        full = run_case(renderer, loss, initial, "contact-thump-resonator", {}, [])
        results = [full]
        for name, fixed, omitted in (
            (
                "contact-thump",
                {"resonance_level": 0},
                [key for key in initial if key.startswith(("resonance_", "tension_"))],
            ),
            (
                "contact-resonator",
                {"thump_level": 0},
                [
                    "thump_pitch_drop_octaves",
                    "thump_pitch_fall_seconds",
                    "thump_decay_seconds",
                ],
            ),
        ):
            results.append(
                run_case(renderer, loss, full["parameters"], name, fixed, omitted)
            )
        summary = dict(
            reference=renderer.metadata,
            start=str(START.relative_to(ROOT)),
            objective=loss.specification,
            results=results,
            listening_approved=False,
        )
        (OUTPUT / "study.json").write_text(json.dumps(summary, indent=2))
        finish_study(renderer, reference, OUTPUT)
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
