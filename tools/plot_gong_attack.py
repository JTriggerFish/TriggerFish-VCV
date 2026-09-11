"""Fixed-level pitch slices and band envelopes for a gong attack trial."""

import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import stft, welch

from triggerfish_percussion.audio_io import read_wav


def figure(directory, baseline):
    reference = read_wav(directory / "reference.wav").mono()
    rate = reference.sample_rate
    sources = [
        ("Reference", reference, "#eebc59"),
        ("Model", read_wav(directory / "candidate.wav").mono(), "#68b5ed"),
    ]
    if baseline:
        sources.append(("Before", read_wav(baseline).mono(), "#bc83c4"))
    fig = make_subplots(
        rows=3,
        cols=2,
        subplot_titles=[
            "Early pitched body: 40–200 ms",
            "Established low body: 200–500 ms",
            "90–180 Hz energy",
            "300–450 Hz energy",
            "500–700 Hz energy",
            "3–12 kHz bloom",
        ],
    )
    bands = [(90, 180), (300, 450), (500, 700), (3000, 12000)]
    for name, audio, colour in sources:
        if audio.sample_rate != rate or len(audio.samples) != len(reference.samples):
            raise ValueError("Expected like-for-like sample rates and lengths")
        for col, (a, b) in enumerate(((0.04, 0.2), (0.2, 0.5)), 1):
            segment = audio.samples[round(a * rate) : round(b * rate)]
            f, p = welch(segment, rate, nperseg=min(8192, len(segment)), nfft=32768)
            selected = (f >= 80) & (f <= 900)
            fig.add_trace(
                go.Scatter(
                    x=f[selected],
                    y=10 * np.log10(np.maximum(p[selected], 1e-20)),
                    name=name,
                    legendgroup=name,
                    showlegend=col == 1,
                    line_color=colour,
                ),
                row=1,
                col=col,
            )
        f, t, z = stft(
            audio.samples, rate, nperseg=4096, noverlap=4096 - round(0.01 * rate)
        )
        for i, (lo, hi) in enumerate(bands):
            energy = np.sum(abs(z[(f >= lo) & (f < hi)]) ** 2, axis=0)
            fig.add_trace(
                go.Scatter(
                    x=t,
                    y=10 * np.log10(np.maximum(energy, 1e-20)),
                    name=name,
                    legendgroup=name,
                    showlegend=False,
                    line_color=colour,
                ),
                row=2 + i // 2,
                col=1 + i % 2,
            )
    fig.update_xaxes(title_text="Hz", range=[80, 900], row=1)
    fig.update_yaxes(title_text="dB/Hz", range=[-75, -25], row=1)
    for row in (2, 3):
        fig.update_xaxes(title_text="seconds", range=[0, 1.5], row=row)
        fig.update_yaxes(title_text="dB power", range=[-70, -10], row=row)
    fig.update_layout(
        template="plotly_dark",
        width=1450,
        height=1100,
        title="Gong attack: absolute levels, same analysis for all signals; no normalization",
    )
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--baseline", type=Path)
    args = p.parse_args()
    (args.directory / "attack.plotly.json").write_text(
        figure(args.directory, args.baseline).to_json(), encoding="utf8"
    )
