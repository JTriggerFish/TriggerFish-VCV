"""Read-only onset/bloom diagnostics with causal bands and an exact saved gesture."""

import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import butter, sosfilt

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

BANDS = [(80, 800), (2500, 5000), (5000, 9000), (9000, 14000)]


def measure(audio, rate):
    """10 ms disjoint RMS windows after causal fourth-order band filtering."""
    hop = round(0.01 * rate)
    frames = len(audio) // hop
    t = (np.arange(frames) + 0.5) * hop / rate
    curves, rows = [], []
    for lo, hi in BANDS:
        signal = sosfilt(
            butter(4, [lo, hi], btype="bandpass", fs=rate, output="sos"), audio
        )
        power = np.mean(signal[: frames * hop].reshape(frames, hop) ** 2, axis=1)
        db = 10 * np.log10(np.maximum(power, 1e-30))
        peak = np.argmax(power[(t >= 0.1) & (t < 2)]) + np.flatnonzero(t >= 0.1)[0]

        def window(start, end):
            return float(
                10 * np.log10(max(1e-30, np.mean(power[(t >= start) & (t < end)])))
            )

        rows.append(
            dict(
                band=[lo, hi],
                first_20ms_db=window(0, 0.02),
                early_20_100ms_db=window(0.02, 0.1),
                bloom_300_1000ms_db=window(0.3, 1),
                peak_db=float(db[peak]),
                peak_seconds=float(t[peak]),
            )
        )
        curves.append(db)
    return t, curves, rows


def run():
    output = Path("build/gong-early-highs")
    output.mkdir(parents=True, exist_ok=True)
    source = json.loads(Path("workbench/web/gong_calibration.fit.json").read_text())
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(r, source)
        p = saved.initial
        cases = {"Current gong": p, "Diffusion off": dict(p, bloom_rate=0)}
        for slope in (-36, -28, -24, -12):
            cases[f"Excitation tilt {slope}"] = dict(p, body_brightness=slope)
        cases["Tilt -24, diffusion off"] = dict(p, body_brightness=-24, bloom_rate=0)
        signals = {"Reference": aligned_reference(r, 3)}
        signals.update(
            {name: saved.render(values, 3) for name, values in cases.items()}
        )
        plot = make_subplots(
            rows=2, cols=2, subplot_titles=[f"{a}–{b} Hz" for a, b in BANDS]
        )
        results = {}
        colours = {
            "Reference": "#eebc59",
            "Current gong": "#68b5ed",
            "Diffusion off": "#999999",
            "Excitation tilt -28": "#c889d4",
        }
        for name, signal in signals.items():
            t, curves, results[name] = measure(signal, r.sample_rate)
            if name not in (
                "Reference",
                "Current gong",
                "Excitation tilt -28",
                "Diffusion off",
            ):
                continue
            for i, curve in enumerate(curves):
                plot.add_trace(
                    go.Scatter(
                        x=t,
                        y=curve,
                        name=name,
                        legendgroup=name,
                        showlegend=i == 0,
                        line_color=colours[name],
                    ),
                    row=i // 2 + 1,
                    col=i % 2 + 1,
                )
        plot.update_xaxes(range=[0, 1.5], title_text="seconds")
        plot.update_yaxes(title_text="Causal band RMS, dBFS")
        plot.update_layout(template="plotly_dark", width=1400, height=900)
        report = dict(
            fit_id=source["id"],
            event=source["controls"]["event"],
            reference=source["reference"],
            parameters=cases,
            results=results,
            note="No preset changes, gain matching, centred STFT or noncausal filtering.",
        )
        (output / "audit.json").write_text(json.dumps(report, indent=2))
        (output / "onset.plotly.json").write_text(plot.to_json())
        print(json.dumps(results, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    run()
