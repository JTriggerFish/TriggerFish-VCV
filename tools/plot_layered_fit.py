"""Show actual reference/candidate band envelopes without changing audio gains."""

import argparse
import json
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.audio_io import read_wav
from audit_gong_bloom_separation import BANDS, measure


def run(args):
    signals = {
        "Reference": args.directory / "reference.wav",
        "Candidate": args.directory / "candidate.wav",
    }
    if args.baseline:
        signals["Previous"] = args.baseline / "candidate.wav"
    plot = make_subplots(
        rows=3, cols=2, subplot_titles=[f"{a}–{b} Hz" for a, b in BANDS]
    )
    report = {}
    for (name, path), colour in zip(signals.items(), ["#eebc59", "#68b5ed", "#9a8bac"]):
        audio = read_wav(path).mono()
        t, db, report[name] = measure(audio.samples, audio.sample_rate)
        for i, band in enumerate(db):
            plot.add_trace(
                go.Scatter(
                    x=t,
                    y=band,
                    name=name,
                    legendgroup=name,
                    showlegend=i == 0,
                    line_color=colour,
                ),
                row=i // 2 + 1,
                col=i % 2 + 1,
            )
    plot.update_xaxes(range=[0, args.seconds], title_text="seconds")
    plot.update_yaxes(title_text="dB band power")
    plot.update_layout(
        template="plotly_dark",
        width=1450,
        height=1100,
        title="Gong: reference-anchored body and upper bloom",
        legend=dict(orientation="h", y=-0.12),
    )
    (args.directory / "layered.plotly.json").write_text(plot.to_json())
    (args.directory / "layered-landmarks.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--seconds", type=float, default=2)
    run(parser.parse_args())
