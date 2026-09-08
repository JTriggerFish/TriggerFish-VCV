"""Offline fixed-scale diagnostic; audition stays in the main workbench."""

import json
import argparse
from pathlib import Path
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.utils import PlotlyJSONEncoder
from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.spectral_difference import spectral_difference


def figure(directory, window=8192, hop=1024, seconds=None):
    ref, candidate = [
        read_wav(directory / name).mono() for name in ("reference.wav", "candidate.wav")
    ]
    if ref.sample_rate != candidate.sample_rate:
        raise ValueError("Sample rates differ")
    reference_samples, candidate_samples = ref.samples, candidate.samples
    if seconds is not None:
        if not np.isfinite(seconds) or seconds <= 0:
            raise ValueError("Crop duration must be positive and finite")
        count = round(seconds * ref.sample_rate)
        if count < 1 or count > min(len(reference_samples), len(candidate_samples)):
            raise ValueError("Crop duration is outside the comparison audio")
        # Crop BEFORE the transform; a later tail must not leak into an attack view.
        reference_samples, candidate_samples = (
            reference_samples[:count],
            candidate_samples[:count],
        )
    data = spectral_difference(
        reference_samples, candidate_samples, ref.sample_rate, window, hop
    )
    selected = (data["frequency"] >= 70) & (data["frequency"] <= 15000)
    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        subplot_titles=[
            "Reference",
            "Model — same scale",
            "Model − reference: red = excess, blue = missing",
        ],
    )
    for row, key in enumerate(("reference", "candidate", "difference"), 1):
        diff = key == "difference"
        fig.add_trace(
            go.Heatmap(
                x=data["time"],
                y=data["frequency"][selected],
                z=data[key][selected],
                colorscale="RdBu_r" if diff else "Inferno",
                zmin=-18 if diff else data["maximum_db"] - 65,
                zmax=18 if diff else data["maximum_db"],
                showscale=diff,
            ),
            row=row,
            col=1,
        )
    ticks = [100, 200, 400, 700, 1000, 2000, 4000, 7000, 10000, 15000]
    fig.update_yaxes(
        type="log",
        range=[np.log10(70), np.log10(15000)],
        tickvals=ticks,
        ticktext=[str(x) for x in ticks],
        title_text="Hz",
    )
    fig.update_xaxes(title_text="seconds", row=3, col=1)
    fig.update_layout(
        width=1450,
        height=1100,
        template="plotly_dark",
        title=f"{directory.name} — {window / ref.sample_rate * 1000:.1f} ms window, "
        f"{hop / ref.sample_rate * 1000:.1f} ms hop",
        margin=dict(l=75, r=70, t=70, b=45),
    )
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--window", type=int, default=8192)
    parser.add_argument("--hop", type=int, default=1024)
    parser.add_argument(
        "--seconds", type=float, help="Crop both signals before the STFT"
    )
    args = parser.parse_args()
    directory = args.directory
    (directory / "difference.plotly.json").write_text(
        json.dumps(
            figure(directory, args.window, args.hop, args.seconds).to_plotly_json(),
            cls=PlotlyJSONEncoder,
        ),
        encoding="utf8",
    )
