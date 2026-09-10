"""Plot low-band motion against the reference; never normalize audition audio."""

import argparse
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.low_mode_beating import LowModeBeating


def figure(directory, baseline=None):
    reference = read_wav(directory / "reference.wav").mono()
    loss = LowModeBeating(reference.samples, reference.sample_rate)
    signals = [("Reference", loss.target, "#f2bc53")]
    paths = [("Candidate", directory / "candidate.wav", "#64b6ec")]
    if baseline is not None:
        paths.append(("Before", baseline, "#c790d9"))
    for label, path, colour in paths:
        audio = read_wav(path).mono()
        if audio.sample_rate != reference.sample_rate:
            raise ValueError("Sample rates differ")
        signals.append((label, loss.analyze(audio.samples), colour))
    titles = [
        title
        for lo, hi in loss.bands
        for title in (f"{lo}–{hi} Hz: relative envelope", "Modulation spectrum")
    ]
    fig = make_subplots(rows=4, cols=2, subplot_titles=titles)
    for label, bands, colour in signals:
        for row, band in enumerate(bands, 1):
            envelope = np.asarray(band["relative_envelope"])
            spectrum = np.asarray(band["spectrum"])
            for col, x, y in (
                (1, 0.5 + np.arange(len(envelope)) / 200, envelope),
                (
                    2,
                    np.linspace(0, 100, len(spectrum)),
                    10 * np.log10(np.maximum(spectrum, 1e-10)),
                ),
            ):
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        name=label,
                        line_color=colour,
                        showlegend=row == 1 and col == 1,
                    ),
                    row=row,
                    col=col,
                )
    fig.update_xaxes(range=[0.5, 4.5], title_text="seconds", col=1)
    fig.update_xaxes(
        type="log",
        range=[np.log10(0.5), np.log10(80)],
        tickvals=[0.5, 1, 3, 8, 20, 40, 80],
        title_text="Hz",
        col=2,
    )
    fig.update_yaxes(range=[-70, 0], title_text="dB/Hz", col=2)
    fig.update_layout(
        template="plotly_dark",
        width=1450,
        height=1100,
        title="Low-band motion — exponential trend removed for this diagnostic only",
    )
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    (args.directory / "motion.plotly.json").write_text(
        figure(args.directory, args.baseline).to_json(), encoding="utf8"
    )
