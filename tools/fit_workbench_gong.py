"""Single-reference gong search using the exact UI/Wasm parameter surface."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.trajectory_fit_loss import TrajectoryLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from workbench_fit_report import write_report

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=ROOT / "build/workbench-wasm/site/gong-review"
    )
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument(
        "--resume", type=Path, default=os.environ.get("TF_GONG_FIT_RESUME")
    )
    parser.add_argument(
        "--audit-only",
        action="store_true",
        default=os.environ.get("TF_GONG_FIT_AUDIT_ONLY") == "1",
    )
    parser.add_argument(
        "--stage",
        choices=("all", "contact", "texture"),
        default=os.environ.get("TF_GONG_FIT_STAGE", "all"),
    )
    parser.add_argument(
        "--overrides", default=os.environ.get("TF_GONG_FIT_OVERRIDES", "{}")
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", ROOT)
    try:
        reference = renderer.reference
        onset = round(
            renderer.metadata["reference"]["cell"].get("onset_seconds", 0)
            * renderer.sample_rate
        )
        reference = reference[onset : onset + 6 * renderer.sample_rate]
        reference = np.pad(
            reference, (0, max(0, 6 * renderer.sample_rate - len(reference)))
        )
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        loss = TrajectoryLoss(reference, renderer.sample_rate)
        search = Search(renderer, loss, args.output)
        if args.resume:
            previous = json.loads(args.resume.read_text())
            search.parameters.update(previous["parameters"])
            search.history = previous["history"]
        overrides = json.loads(args.overrides)
        if overrides:
            search.parameters.update(overrides)
            search.history.append(
                dict(
                    stage="explicit review correction",
                    overrides=overrides,
                    parameters=dict(search.parameters),
                )
            )
        if args.audit_only:
            # Round-trip via the UI's actual import validation and adapter.
            fit = renderer.request(command="snapshot", parameters=search.parameters)[
                "fit"
            ]
            expected = search.audio(search.parameters)
            actual = renderer.decode(
                renderer.request(command="renderSnapshot", fit=fit, seconds=6)["pcm"]
            )
            if not np.array_equal(actual, expected):
                raise RuntimeError(
                    "Saved UI fit does not reproduce the optimized audio"
                )
            search.save()
            heldout = [
                loss.diagnostics(
                    renderer.render(
                        search.parameters,
                        6,
                        renderer.metadata["event"]["seed"] + offset,
                    )
                )
                for offset in (1, 2, 3)
            ]
            (args.output / "heldout-seeds.json").write_text(
                json.dumps(heldout, indent=2)
            )
            write_report(args.output, renderer=renderer)
            print("Saved fit round-trip: sample-identical", flush=True)
            return
        baseline = search.audio(search.parameters)
        if not args.resume or not (args.output / "baseline.wav").exists():
            write_wav(
                args.output / "baseline.wav",
                AudioBuffer(baseline, renderer.sample_rate),
            )
        print(json.dumps(dict(baseline=loss.diagnostics(baseline))), flush=True)
        search.save()
        if args.stage == "texture":
            search.explore_texture()
            write_report(args.output, renderer=renderer)
            return
        levels = {f"resolved_level_{index}": (-30, 24) for index in range(17)}
        if args.stage == "all":
            search.stage("observed spectrum", levels, args.iterations)
        dynamics = {
            "bloom_rate": (0.1, 8),
            "body_brightness": (-72, -6),
            "body_excitation_centre": (300, 3000),
            "field_turbulence": (0.05, 1),
            "field_turbulence_slope": (0.05, 1),
            "field_turbulence_centre": (600, 3500),
            "field_packet_spread": (0.1, 8),
            "field_phase_bandwidth": (0, 2),
            "field_exchange": (0, 1),
            "bloom_phase_diffusion": (0, 1),
        }
        if args.stage == "all":
            search.stage("transport and texture", dynamics, args.iterations)
            search.stage(
                "two-point intrinsic loss",
                {
                    "body_decay_seconds_0": (2, 20),
                    "body_decay_seconds_7": (0.1, 12),
                    "bloom_rate": (0.1, 8),
                },
                args.iterations,
            )
            search.stage("observation refinement", levels, args.iterations)
        search.stage(
            "source and radiation",
            {
                "direct_gain": (0, 0.04),
                "impact_tone_noise": (0, 0.8),
                "impact_width": (0.5, 4),
                "impact_chirp_pitch": (0.05, 0.4),
                "body_high_cut": (6000, 22000),
                "direct_high_cut": (1000, 6000),
            },
            args.iterations,
        )
        # Hold out random realizations; don't optimize each seed's phases.
        seeds = [renderer.metadata["event"]["seed"] + offset for offset in (1, 2, 3)]
        diagnostics = [
            loss.diagnostics(renderer.render(search.parameters, 6, seed))
            for seed in seeds
        ]
        (args.output / "heldout-seeds.json").write_text(
            json.dumps(diagnostics, indent=2)
        )
        write_report(args.output, renderer=renderer)
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
