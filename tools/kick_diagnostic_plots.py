"""Offline source/shape plots; no extra audition server or public report page."""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.transforms import StftConfig, stft
from triggerfish_percussion.band_region_audit import BandRegionAudit, BANDS, REGIONS

COLOURS = dict(
    reference="#fcbf49",
    published="#06b6d4",
    thump_level="#93c5fd",
    resonance_level="#f472b6",
    **{
        "without-source-noise": "#a3a3a3",
        "coverage-bypass": "#4ade80",
        "coverage-lowpass": "#c084fc",
        "measured-lowpass": "#fb7185",
    },
)


def plot_diagnosis(directory):
    reference = read_wav(directory / "reference.wav").mono()
    names = [
        "reference",
        "published",
        "thump_level",
        "resonance_level",
        "without-source-noise",
    ]
    names += [
        name
        for name in ("coverage-bypass", "coverage-lowpass", "measured-lowpass")
        if (directory / f"{name}.wav").exists()
    ]
    figure = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "Attack: 0–30 ms · source contributions",
            "Early decay: 30–100 ms",
            "Band/time error of published preset",
            "Diagnostic trials · 0–30 ms",
        ),
    )
    for name in names:
        audio = read_wav(directory / f"{name}.wav").mono()
        spectrum = stft(audio.samples, audio.sample_rate, StftConfig(2048, 128))
        power = spectrum.power
        for col, (a, b) in enumerate(((0, 0.03), (0.03, 0.1)), 1):
            selected = (spectrum.times_seconds >= a) & (spectrum.times_seconds < b)
            y = 10 * np.log10(np.maximum(power[:, selected].mean(axis=1), 1e-14))
            if name in (
                "reference",
                "published",
                "thump_level",
                "resonance_level",
                "without-source-noise",
            ):
                figure.add_trace(
                    go.Scatter(
                        x=spectrum.frequencies_hz,
                        y=y,
                        name=name,
                        showlegend=col == 1,
                        line=dict(color=COLOURS[name]),
                    ),
                    row=1,
                    col=col,
                )
            if col == 1 and name in (
                "reference",
                "published",
                "coverage-bypass",
                "coverage-lowpass",
                "measured-lowpass",
            ):
                figure.add_trace(
                    go.Scatter(
                        x=spectrum.frequencies_hz,
                        y=y,
                        name=name,
                        line=dict(color=COLOURS[name]),
                        showlegend=name not in ("reference", "published"),
                    ),
                    row=2,
                    col=2,
                )
    audit = BandRegionAudit(reference.samples, reference.sample_rate)
    measured = audit.measure(read_wav(directory / "published.wav").mono().samples)
    errors = np.array(measured["error_db"])
    errors[~np.array(measured["evaluated"])] = np.nan
    figure.add_trace(
        go.Heatmap(
            z=errors.tolist(),
            zmin=-15,
            zmax=15,
            colorscale="RdBu",
            reversescale=True,
            x=[f"{a*1000:g}–{b*1000:g}ms" for a, b in REGIONS],
            y=[f"{a}–{b}Hz" for a, b in BANDS],
            showscale=True,
            colorbar=dict(title="Δ dB", x=0.45, y=0.2, len=0.35, thickness=12),
        ),
        row=2,
        col=1,
    )
    for row, col in ((1, 1), (1, 2), (2, 2)):
        figure.update_xaxes(
            type="log",
            range=[np.log10(100), np.log10(10000)],
            title="Hz",
            row=row,
            col=col,
        )
        figure.update_yaxes(
            range=[-115, -20], title="Spectral power dBFS", row=row, col=col
        )
    figure.update_layout(
        height=980,
        width=1500,
        template="plotly_dark",
        title="Kick diagnosis · fixed gain · trials are not published",
    )
    figure.write_html(directory / "analysis.html", include_plotlyjs=True)
