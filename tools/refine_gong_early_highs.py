"""Fit only initial excitation tilt and the existing high T60 endpoint.

Fixed causal-band time regions keep the quiet onset and bloom independently
represented. No modal placement, observation gains, EQ or level matching.
"""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from audit_gong_early_highs import measure
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def validate(saved, reference, best, output):
    """Inspect the full tail and an unoptimized seed, not just the search score."""
    plot = make_subplots(
        rows=2, cols=2, subplot_titles=["80–800 Hz", "2.5–5 kHz", "5–9 kHz", "9–14 kHz"]
    )
    signals = {
        "Reference": reference,
        "Original": saved.render(saved.initial, 8),
        "Revised": saved.render(best, 8),
        "Revised, holdout seed": saved.render(best, 8, 73519),
    }
    colours = ["#eebc59", "#68b5ed", "#c889d4", "#80c6a4"]
    rows = {}
    for (name, signal), colour in zip(signals.items(), colours):
        t, curves, rows[name] = measure(signal, saved.sample_rate)
        for i, curve in enumerate(curves):
            plot.add_trace(
                go.Scatter(
                    x=t,
                    y=curve,
                    name=name,
                    legendgroup=name,
                    showlegend=i == 0,
                    line_color=colour,
                ),
                row=i // 2 + 1,
                col=i % 2 + 1,
            )
    plot.update_xaxes(range=[0, 5], title_text="seconds")
    plot.update_yaxes(range=[-125, -10], title_text="Causal band RMS, dBFS")
    plot.update_layout(template="plotly_dark", width=1400, height=900)
    (output / "validation.plotly.json").write_text(plot.to_json())
    return rows


def regions(audio, rate):
    t, curves, _ = measure(audio, rate)
    power = 10 ** (np.asarray(curves) / 10)
    windows = [(0.02, 0.1), (0.1, 0.3), (0.3, 1), (1, 2), (2, 3)]
    return np.asarray(
        [
            10 * np.log10(np.maximum(1e-30, power[:, (t >= a) & (t < b)].mean(axis=1)))
            for a, b in windows
        ]
    ).T


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fit", type=Path, default=Path("workbench/web/gong_calibration.fit.json")
    )
    args = parser.parse_args()
    output = Path("build/gong-early-refinement")
    output.mkdir(parents=True, exist_ok=True)
    fit = json.loads(args.fit.read_text(encoding="utf-8"))
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(renderer, fit)
        reference = aligned_reference(renderer, 3)
        target = regions(reference, renderer.sample_rate)
        seeds = (fit["controls"]["event"]["seed"], 1982)
        cases = [("Original", saved.initial)]
        cases.extend(
            (
                f"Tilt {tilt}, upper T60 {decay}",
                dict(saved.initial, body_brightness=tilt, body_decay_seconds_7=decay),
            )
            for tilt, decay in product((-30, -28, -26), (1.25, 1.5, 1.75))
        )
        rows = []
        for name, parameters in cases:
            observations = [
                regions(saved.render(parameters, 3, seed), renderer.sample_rate)
                for seed in seeds
            ]
            errors = np.asarray(observations) - target
            # Fit low body + the two high bands; the user's scooped upper mids
            # are recorded but not corrected by distorting unrelated controls.
            score = float(np.sqrt(np.mean(errors[:, [0, 2, 3], :] ** 2)))
            rows.append(
                dict(
                    name=name,
                    parameters=parameters,
                    score=score,
                    regions=np.mean(observations, axis=0).tolist(),
                    errors=errors.tolist(),
                )
            )
            print(json.dumps(dict(name=name, score=score)), flush=True)
        best = min(rows, key=lambda row: row["score"])
        candidate = saved.snapshot(best["parameters"], "Gong — quiet onset and bloom")
        report = dict(
            source=fit,
            target=target.tolist(),
            seeds=seeds,
            rows=rows,
            selected=best["name"],
            renderer_sha256=renderer.metadata["rendererSha256"],
        )
        report["validation"] = validate(
            saved, aligned_reference(renderer, 8), best["parameters"], output
        )
        (output / "candidate.fit.json").write_text(json.dumps(candidate, indent=2))
        (output / "report.json").write_text(json.dumps(report, indent=2))
        print(
            json.dumps(
                dict(
                    selected=best["name"],
                    target=target.tolist(),
                    regions=best["regions"],
                )
            )
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    run()
