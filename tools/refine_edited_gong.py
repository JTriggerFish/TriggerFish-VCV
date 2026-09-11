"""Refine a saved gong with shared dynamics/texture only; never fit a ridge.

Reference gain, saved gesture, all modal bars and frequencies stay fixed.
The hybrid target explicitly retains the edited low body and compares the
upper sound to the recording. Trials are evaluated at two separate seeds.
"""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from audit_gong_bloom_separation import BANDS, measure
from gong_texture_comparison import GongTextureComparison


def plot_signals(signals, rate, output):
    fig = make_subplots(
        rows=3, cols=2, subplot_titles=[f"{a}–{b} Hz" for a, b in BANDS]
    )
    rows = {}
    for (name, audio), colour in zip(
        signals.items(), ("#eebc59", "#bc83c4", "#68b5ed")
    ):
        t, db, rows[name] = measure(audio, rate)
        for i, values in enumerate(db):
            fig.add_trace(
                go.Scatter(
                    x=t.tolist(),
                    y=values.tolist(),
                    name=name,
                    line_color=colour,
                    legendgroup=name,
                    showlegend=i == 0,
                ),
                row=i // 2 + 1,
                col=i % 2 + 1,
            )
    fig.update_xaxes(range=[0, 6], title_text="seconds")
    fig.update_yaxes(title_text="dB band power")
    fig.update_layout(template="plotly_dark", width=1450, height=1050)
    (output / "envelopes.plotly.json").write_text(fig.to_json())
    (output / "landmarks.json").write_text(json.dumps(rows, indent=2))


class Comparison:
    def __init__(self, reference, body, rate):
        self.spectral = GongTextureComparison(reference, body, rate)
        self.temporal = SpectralBloomLoss(reference, rate)
        centres = np.sqrt(self.temporal.edges[:-1] * self.temporal.edges[1:])
        mix = np.clip(np.log2(centres / 800) / np.log2(1800 / 800), 0, 1)[:, None]
        mix = mix * mix * (3 - 2 * mix)
        self.target = (1 - mix) * self.temporal.db(
            self.temporal.power(body)
        ) + mix * self.temporal.target
        self.low = centres < 800

    def metrics(self, audio):
        error = self.temporal.db(self.temporal.power(audio)) - self.target
        # Separate a constant level mismatch from envelope shape. No audition
        # normalization: the incoming bars/gains are fixed and all raw errors
        # remain reported. Do not shorten T60 to turn down a loud upper packet.
        bias = error.mean(axis=1, keepdims=True)
        timing = float(np.sqrt(np.mean((error - bias) ** 2)))
        body = float(np.sqrt(np.mean(error[self.low] ** 2)))
        spectral = self.spectral.score(audio)
        return dict(
            score=float(np.sqrt(timing**2 + body**2)),
            spectral=spectral,
            timing=timing,
            raw_envelope_db=float(np.sqrt(np.mean(error**2))),
            body_envelope_db=body,
            **self.spectral.measure(audio),
        )


def run(args):
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    source = json.loads(args.fit.read_text())
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(r, source)
        seeds = (source["controls"]["event"]["seed"], 1982)
        reference = aligned_reference(r, 6)
        before = [saved.render(saved.initial, 6, seed) for seed in seeds]
        targets = [Comparison(reference, audio, r.sample_rate) for audio in before]
        p = dict(saved.initial)
        if p["output_eq_enabled"] or p["field_phase_bandwidth"]:
            raise ValueError("Expected the user's EQ-free movement-only starting point")
        history, cache = [], {}

        def evaluate(parameters):
            key = tuple(parameters.items())
            if key not in cache:
                metrics = [
                    target.metrics(saved.render(parameters, 6, seed))
                    for target, seed in zip(targets, seeds)
                ]
                cache[key] = dict(
                    score=float(np.mean([m["score"] for m in metrics])), metrics=metrics
                )
            return cache[key]

        baseline = evaluate(p)
        print(json.dumps(dict(baseline=baseline)), flush=True)
        plot_signals(
            {"Reference": reference, "Your edit": before[0]}, r.sample_rate, output
        )
        write_wav(output / "reference.wav", AudioBuffer(reference, r.sample_rate))
        write_wav(output / "edited.wav", AudioBuffer(before[0], r.sample_rate))
        if args.inspect:
            return
        stages = [
            (
                "shared dynamics",
                {
                    "body_decay_seconds_7": (1.3, 1.7, 2.14, 2.7, 3.5),
                    "body_decay_seconds_0": (7, 8.5, 10, 12, 14),
                    "bloom_rate": (2.5, 3.2, 4, 5, 6.5),
                    "bloom_energy_acceleration": (0.015, 0.03, 0.05, 0.085, 0.14),
                },
            ),
            (
                "shared texture",
                {
                    "field_motion_depth": (0.75, 1, 1.25, 1.5, 1.8),
                    "field_motion_rate": (60, 100, 150, 200),
                    "field_packet_spread": (1, 1.4, 1.8, 2.3),
                },
            ),
        ]
        if args.fine:
            # Refine the shared dynamics jointly, retaining the edited texture.
            # The coarse movement-depth change had negligible timing benefit
            # and moved ridge contrast away from the recording.
            rows = [dict(parameters=p, **evaluate(p))]
            if args.reuse:
                previous = json.loads(args.reuse.read_text())
                rows = previous["history"][0]["trials"]
                if (
                    previous["source_id"] != source["id"]
                    or previous["event"] != source["controls"]["event"]
                    or previous["renderer_sha256"] != r.metadata["rendererSha256"]
                    or rows[0]["parameters"] != p
                ):
                    raise ValueError(
                        "Cached scan does not describe this snapshot/renderer"
                    )
            grid = (
                []
                if args.reuse
                else product((2.4, 2.8, 3.2, 3.6, 4), (1.5, 1.75, 2, 2.25))
            )
            for rate, decay in grid:
                trial = dict(p, bloom_rate=rate, body_decay_seconds_7=decay)
                row = dict(parameters=trial, **evaluate(trial))
                rows.append(row)
                print(
                    json.dumps(dict(rate=rate, decay=decay, score=row["score"])),
                    flush=True,
                )
            # Do not improve a pooled timing score by worsening fine-band decay.
            # The unconstrained winner visibly over-shortened the highest tail.
            eligible = [
                row
                for row in rows
                if all(
                    m["upper_shape_db"] <= b["upper_shape_db"] + 1e-8
                    for m, b in zip(row["metrics"], baseline["metrics"])
                )
            ]
            winner = min(eligible, key=lambda row: row["score"])
            p = winner["parameters"]
            history.append(
                dict(
                    stage="joint shared dynamics; texture retained",
                    trials=rows,
                    acceptance="no worse fine-band upper decay-shape RMS on either seed",
                )
            )
        for stage, grid in ([] if args.fine else stages):
            for key, values in grid.items():
                rows = [(p[key], evaluate(p))]
                for value in values:
                    candidate = dict(p, **{key: value})
                    row = evaluate(candidate)
                    rows.append((value, row))
                    print(
                        json.dumps(dict(stage=stage, key=key, value=value, **row)),
                        flush=True,
                    )
                value, metrics = min(rows, key=lambda row: row[1]["score"])
                history.append(
                    dict(stage=stage, key=key, before=p[key], after=value, trials=rows)
                )
                p[key] = value
                (output / "progress.json").write_text(
                    json.dumps(dict(parameters=p, history=history), indent=2)
                )
        after = saved.render(p, 6)
        candidate = saved.snapshot(p, "Gong — edited series, refined dynamics")
        # Exact snapshot replay, not just a matching set of printed parameters.
        replay = r.decode(
            r.request(command="renderSnapshot", fit=candidate, seconds=6)["pcm"]
        )
        if not np.array_equal(after, replay):
            raise ValueError("Exported candidate does not reproduce")
        if any(
            p[k] != v for k, v in saved.initial.items() if k.startswith("resolved_")
        ):
            raise ValueError("The modal series was changed")
        (output / "candidate.fit.json").write_text(json.dumps(candidate, indent=2))
        (output / "original.fit.json").write_text(json.dumps(source, indent=2))
        write_wav(output / "candidate.wav", AudioBuffer(after, r.sample_rate))
        plot_signals(
            {"Reference": reference, "Your edit": before[0], "Candidate": after},
            r.sample_rate,
            output,
        )
        report = dict(
            source_id=source["id"],
            event=source["controls"]["event"],
            seeds=seeds,
            renderer_sha256=r.metadata["rendererSha256"],
            target=targets[0].spectral.specification,
            objective="hypot(equal-band equal-time envelope shape RMS, preserved low-body envelope RMS); average two seeds. Constant band-level error reported separately, not fitted through damping",
            fixed="all individual mode parameters, EQ off, phase blur off, reference gain and gesture",
            baseline=baseline,
            final=evaluate(p),
            history=history,
            changed={
                k: [saved.initial[k], v] for k, v in p.items() if saved.initial[k] != v
            },
        )
        (output / "audit.json").write_text(json.dumps(report, indent=2))
        print(
            json.dumps(dict(changed=report["changed"], final=report["final"])),
            flush=True,
        )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inspect", action="store_true")
    parser.add_argument("--fine", action="store_true")
    parser.add_argument("--reuse", type=Path)
    run(parser.parse_args())
