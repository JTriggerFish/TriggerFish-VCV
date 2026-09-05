"""Kick-specific displays: retain attack timing separately from bass resolution."""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d

from triggerfish_percussion.transforms import StftConfig, stft


def detail_plot(reference, candidate):
    figure = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            "Waveform — first 180 ms, identical axes",
            "12 ms RMS envelope",
            "Attack spectrum — 512-sample window",
            "Low spectrum — 8192-sample window",
        ],
    )
    for name, audio, colour, dash in (
        ("Reference", reference, "#fcbf49", "solid"),
        ("Candidate", candidate, "#06d6a0", "dash"),
    ):
        rate, samples = audio.sample_rate, audio.samples
        time = np.arange(len(samples)) / rate
        limit = round(0.18 * rate)
        figure.add_trace(
            go.Scatter(
                x=time[:limit:8],
                y=samples[:limit:8],
                name=name,
                line=dict(color=colour, dash=dash),
            ),
            row=1,
            col=1,
        )
        power = gaussian_filter1d(samples**2, 0.012 * rate)
        figure.add_trace(
            go.Scatter(
                x=time[::88],
                y=10 * np.log10(np.maximum(power[::88], 1e-12)),
                name=name,
                showlegend=False,
                line=dict(color=colour, dash=dash),
            ),
            row=1,
            col=2,
        )
        for column, size in ((1, 512), (2, 8192)):
            value = stft(samples, rate, StftConfig(size, 128))
            for start, end, opacity in (
                (0, 0.03, 1),
                (0.03, 0.1, 0.65),
                (0.1, 0.25, 0.35),
            ):
                selected = (value.times_seconds >= start) & (value.times_seconds < end)
                spectrum = value.power[:, selected].mean(axis=1)
                figure.add_trace(
                    go.Scatter(
                        x=value.frequencies_hz,
                        y=10 * np.log10(np.maximum(spectrum, 1e-14)),
                        name=f"{name} {start*1000:g}–{end*1000:g} ms",
                        opacity=opacity,
                        showlegend=column == 2,
                        line=dict(color=colour, dash=dash),
                    ),
                    row=2,
                    col=column,
                )
    figure.update_xaxes(title="Seconds", row=1)
    figure.update_yaxes(title="Amplitude", row=1, col=1)
    figure.update_yaxes(title="dBFS", range=[-100, 0], row=1, col=2)
    figure.update_xaxes(
        type="log", range=[np.log10(250), np.log10(16000)], title="Hz", row=2, col=1
    )
    figure.update_xaxes(range=[20, 300], title="Hz", row=2, col=2)
    figure.update_yaxes(title="Spectral power (dBFS)", range=[-100, -10], row=2)
    figure.update_layout(
        height=700, template="plotly_dark", margin=dict(l=60, r=30, t=45, b=40)
    )
    return figure.to_html(full_html=False, include_plotlyjs=False)
