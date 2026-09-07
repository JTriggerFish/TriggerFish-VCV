"""Offline shared-scale plots: each trial is compared to the real reference."""

import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.transforms import StftConfig, stft


def draw_comparison(directory):
    rows = json.loads((directory / "trials.json").read_text(encoding="utf8"))
    reference = read_wav(directory / "reference.wav").mono()
    winners = []
    for objective in ("region", "mel", "mel_a", "jtfs"):
        trials = [r for r in rows if r["name"].rsplit("-direct-", 1)[0] == objective]
        if trials:
            winners.append(min(trials, key=lambda r: r["fit"]["after"]))
    config = StftConfig(2048, 128)
    target = stft(reference.samples, reference.sample_rate, config)
    ceiling = float(10 * np.log10(max(target.power.max(), 1e-20)))
    figure = make_subplots(
        rows=len(winners),
        cols=2,
        subplot_titles=[
            title
            for r in winners
            for title in ("Reference", r["name"] + " — not approved")
        ],
        vertical_spacing=0.07,
        horizontal_spacing=0.05,
    )
    for i, row in enumerate(winners, 1):
        candidate = read_wav(directory / row["name"] / "candidate.wav").mono()
        for j, audio in enumerate((reference, candidate), 1):
            value = stft(audio.samples, audio.sample_rate, config)
            bins = (value.frequencies_hz >= 20) & (value.frequencies_hz <= 8000)
            times = value.times_seconds <= 0.65
            power = value.power[np.ix_(bins, times)]
            figure.add_trace(
                go.Heatmap(
                    x=value.times_seconds[times],
                    y=value.frequencies_hz[bins],
                    z=10 * np.log10(np.maximum(power, 1e-20)),
                    coloraxis="coloraxis",
                ),
                i,
                j,
            )
    figure.update_xaxes(title_text="seconds", range=[0, 0.65])
    figure.update_yaxes(
        type="log",
        range=[np.log10(20), np.log10(8000)],
        tickvals=[30, 60, 120, 250, 500, 1000, 2000, 4000, 8000],
    )
    figure.update_layout(
        template="plotly_dark",
        height=1100,
        width=1500,
        title="Matched loss trials: fixed gain and reference-only colour scale",
        coloraxis=dict(
            colorscale="Magma",
            cmin=ceiling - 70,
            cmax=ceiling,
            colorbar=dict(title="dBFS"),
        ),
    )
    figure.write_html(directory / "comparison.html", include_plotlyjs=True)


if __name__ == "__main__":
    import sys

    draw_comparison(Path(sys.argv[1]))
