"""Separate pitched-body motion from a migrating spectral-energy envelope.

Read-only ablations at the saved reference gesture; no preset edits, independent
audio normalization, synthetic delay, or fitting. Timing thresholds are relative
to each band's own peak for diagnosis, never applied to audition audio.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d
from scipy.signal import stft

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

BANDS = [(80, 250), (300, 800), (800, 2500), (2500, 5000), (5000, 9000), (9000, 14000)]


def measure(audio, rate):
    f, t, z = stft(audio, rate, nperseg=4096, noverlap=4096 - round(0.005 * rate))
    power = abs(z) ** 2
    rows, envelopes = [], []
    for lo, hi in BANDS:
        p = gaussian_filter1d(power[(f >= lo) & (f < hi)].sum(axis=0), 4)
        region = (t >= 0) & (t < 2)
        peak = int(np.flatnonzero(region)[np.argmax(p[region])])
        peak_power = p[peak]
        crossings = np.flatnonzero((p >= peak_power * 0.5) & (t <= t[peak]))
        rows.append(
            dict(
                band_hz=[lo, hi],
                peak_seconds=float(t[peak]),
                first_half_power_seconds=float(t[crossings[0]]),
                peak_db=float(10 * np.log10(max(peak_power, 1e-20))),
            )
        )
        envelopes.append(10 * np.log10(np.maximum(p, 1e-12)))
    return t, np.array(envelopes), rows


def run(args):
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        base = dict(r.initial)
        cases = {
            "Current": base,
            "Movement off": dict(base, field_motion_depth=0),
            "Both motion and wander off": dict(
                base, field_motion_depth=0, field_wander_hz=0
            ),
            "Diffusion off": dict(base, bloom_rate=0),
            "Faster diffusion (diagnostic only)": dict(base, bloom_rate=12),
        }
        signals = {"Reference": aligned_reference(r, 6)}
        signals.update({name: r.render(p, 6) for name, p in cases.items()})
        plot = make_subplots(
            rows=3, cols=2, subplot_titles=[f"{a}–{b} Hz" for a, b in BANDS]
        )
        report = dict(
            parameters=base,
            event=r.metadata["event"],
            reference=r.metadata["reference"],
            fft_size=4096,
            hop_seconds=0.005,
            power_smoothing_sigma_seconds=0.02,
            results={},
        )
        colours = ["#eebc59", "#68b5ed", "#ed8ac5", "#89c988", "#a6a6a6", "#bf96fc"]
        for (name, audio), colour in zip(signals.items(), colours):
            t, db, rows = measure(audio, r.sample_rate)
            report["results"][name] = rows
            for i, y in enumerate(db):
                plot.add_trace(
                    go.Scatter(
                        x=t,
                        y=y,
                        name=name,
                        legendgroup=name,
                        showlegend=i == 0,
                        line_color=colour,
                    ),
                    row=i // 2 + 1,
                    col=i % 2 + 1,
                )
            write_wav(
                args.output / (name.split(" (")[0].lower().replace(" ", "-") + ".wav"),
                AudioBuffer(audio, r.sample_rate),
            )
        plot.update_xaxes(range=[0, 2], title_text="seconds")
        plot.update_yaxes(title_text="dB band power")
        plot.update_layout(
            template="plotly_dark",
            width=1450,
            height=1100,
            title="Gong: stable body plus late sizzle, or continuous spectral travel?",
            legend=dict(orientation="h", y=-0.12),
        )
        (args.output / "separation.plotly.json").write_text(plot.to_json())
        (args.output / "audit.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report["results"], indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
