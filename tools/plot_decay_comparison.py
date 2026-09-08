"""Offline absolute band-envelope comparison; no candidate level normalization."""

import argparse
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.utils import PlotlyJSONEncoder
from scipy.ndimage import gaussian_filter1d

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss


def figure(directory):
    """Show matching scales and a reference-derived floor in every band."""
    reference, candidate = [
        read_wav(directory / name).mono() for name in ("reference.wav", "candidate.wav")
    ]
    if reference.sample_rate != candidate.sample_rate:
        raise ValueError("Sample rates differ")
    if reference.samples.shape != candidate.samples.shape:
        raise ValueError("Comparison durations differ")
    loss = BandDecayShapeLoss(reference.samples, reference.sample_rate)
    powers = [loss.power(audio.samples) for audio in (reference, candidate)]
    titles = [f"{low}–{high} Hz" for low, high in loss.bands]
    result = make_subplots(rows=3, cols=2, subplot_titles=titles, shared_xaxes=True)
    sigma = 0.02 * reference.sample_rate / loss.config.hop_samples
    for index, (ref, synth) in enumerate(zip(*powers)):
        row, col = index // 2 + 1, index % 2 + 1
        for name, power, colour in (
            ("Reference", ref, "#efb94d"),
            ("TriggerFish", synth, "#68c4f2"),
        ):
            # Smooth power, not dB; identical 20 ms sigma for both recordings.
            db = 10 * np.log10(np.maximum(gaussian_filter1d(power, sigma), loss.floor))
            result.add_trace(
                go.Scatter(
                    x=loss.times,
                    y=db,
                    name=name,
                    legendgroup=name,
                    showlegend=index == 0,
                    line=dict(color=colour),
                ),
                row=row,
                col=col,
            )
    result.update_xaxes(
        title_text="Seconds", range=[0, loss.frames / loss.rate], showticklabels=True
    )
    result.update_yaxes(title_text="Band power (dB)")
    result.update_layout(
        title="Absolute decay envelopes — unchanged audio levels",
        template="plotly_dark",
        width=1450,
        height=1000,
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    payload = figure(args.directory).to_plotly_json()
    (args.directory / "decay.plotly.json").write_text(
        json.dumps(payload, cls=PlotlyJSONEncoder), encoding="utf8"
    )
