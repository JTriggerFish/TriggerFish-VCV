"""Private offline inspection image; the user's audition remains the workbench."""

import json
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.utils import PlotlyJSONEncoder
from scipy.signal import butter, sosfilt

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.power_envelope import smoothed_power
from triggerfish_percussion.transforms import StftConfig, stft


def region_spectrum(samples, rate, start, end):
    """Crop before windowing: later energy must not leak into attack plots."""
    segment = samples[round(start * rate) : round(end * rate)]
    if not len(segment):
        raise ValueError("Spectrum region is outside the audio")
    value = stft(segment, rate, StftConfig(4096, 512))
    selected = value.times_seconds < len(segment) / rate
    return value.frequencies_hz, value.power[:, selected].mean(axis=1)


def figure(directory):
    reference_path = directory / "reference.wav"
    if not reference_path.exists():
        reference_path = directory.parent / "reference.wav"
    audio = [
        read_wav(path).mono() for path in (reference_path, directory / "candidate.wav")
    ]
    rate = audio[0].sample_rate
    if audio[1].sample_rate != rate or len(audio[0].samples) != len(audio[1].samples):
        raise ValueError("Comparison audio must share a sample rate and duration")
    long = [stft(x.samples, rate, StftConfig(4096, 512)) for x in audio]
    maximum = 10 * np.log10(max(long[0].power.max(), 1e-20))
    fig = make_subplots(
        rows=3,
        cols=2,
        subplot_titles=[
            "Reference — fixed scale",
            "Candidate — same scale",
            "Initial 120 ms spectrum",
            "0.12–0.5 s spectrum",
            "Bass / low-mid power envelopes",
            "Mid / upper power envelopes",
        ],
    )
    for column, value in enumerate(long, 1):
        chosen = (value.frequencies_hz >= 40) & (value.frequencies_hz <= 16000)
        fig.add_trace(
            go.Heatmap(
                x=value.times_seconds.tolist(),
                y=value.frequencies_hz[chosen].tolist(),
                z=(10 * np.log10(np.maximum(value.power[chosen], 1e-20))).tolist(),
                zmin=maximum - 65,
                zmax=maximum,
                colorscale="Inferno",
                showscale=False,
            ),
            row=1,
            col=column,
        )
    for j, (label, colour) in enumerate(
        (("Reference", "#fcbf49"), ("Candidate", "#06d6a0"))
    ):
        for column, (start, end) in enumerate(((0, 0.12), (0.12, 0.5)), 1):
            frequencies, power = region_spectrum(audio[j].samples, rate, start, end)
            fig.add_trace(
                go.Scatter(
                    x=frequencies.tolist(),
                    y=(10 * np.log10(np.maximum(power, 1e-20))).tolist(),
                    name=label,
                    legendgroup=label,
                    showlegend=column == 1,
                    line=dict(color=colour, width=1),
                ),
                row=2,
                col=column,
            )
        for band_index, (low, high) in enumerate(
            ((40, 250), (250, 1000), (1000, 4000), (4000, 16000))
        ):
            sos = butter(
                2,
                [low, min(high, rate * 0.49)],
                btype="bandpass",
                fs=rate,
                output="sos",
            )
            power = smoothed_power(sosfilt(sos, audio[j].samples), 0.012 * rate)[::256]
            fig.add_trace(
                go.Scatter(
                    x=(np.arange(len(power)) * 256 / rate).tolist(),
                    y=(10 * np.log10(np.maximum(power, 1e-12))).tolist(),
                    name=f"{label} {low}–{high} Hz",
                    line=dict(
                        color=colour, dash="solid" if band_index % 2 == 0 else "dot"
                    ),
                ),
                row=3,
                col=1 + band_index // 2,
            )
    fig.update_yaxes(type="log", range=[np.log10(40), np.log10(16000)], row=1)
    fig.update_xaxes(type="log", range=[np.log10(40), np.log10(16000)], row=2)
    fig.update_yaxes(range=[maximum - 75, maximum + 5], row=2)
    fig.update_yaxes(range=[-100, 0], row=3)
    fig.update_layout(
        width=1450,
        height=1100,
        template="plotly_dark",
        title=directory.name,
        margin=dict(l=55, r=20, t=65, b=30),
        legend=dict(orientation="h", y=-0.07),
    )
    return fig


if __name__ == "__main__":
    directory = Path(sys.argv[1])
    (directory / "inspection.plotly.json").write_text(
        json.dumps(figure(directory).to_plotly_json(), cls=PlotlyJSONEncoder),
        encoding="utf8",
    )
