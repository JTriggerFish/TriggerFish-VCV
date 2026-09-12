"""Fixed-level crash comparisons: separate spectral, decay and texture errors."""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import stft, welch
from scipy.ndimage import gaussian_filter1d

from crash_beating_common import CrashObjective
from triggerfish_percussion.band_decay_shape_loss import BandDecayShapeLoss


def ridge_contrast(audio, rate):
    """Descriptive line-versus-wash statistic after removing broad spectral colour.

    One-second Hann spectra, locally whitened by a 50-Hz Gaussian power average.
    Negative dB means more concentrated ridges; not an acceptance loss.
    """
    rows = []
    for start in (0.1, 0.6, 1.2):
        segment = audio[round(start * rate) : round((start + 1) * rate)]
        f, power = welch(segment, rate, nperseg=len(segment), noverlap=0)
        floor = max(float(power.max()) * 1e-6, 1e-25)
        shape = gaussian_filter1d(power, 50 / (f[1] - f[0]))
        whitened = np.maximum(power, floor) / np.maximum(shape, floor)
        row = []
        for low, high in ((1500, 3000), (3000, 6000), (6000, 12000)):
            values = whitened[(f >= low) & (f < high) & (shape > floor * 10)]
            row.append(
                float(10 * np.log10(np.exp(np.mean(np.log(values))) / np.mean(values)))
                if len(values)
                else None
            )
        rows.append(row)
    return dict(
        regions_seconds=[[0.1, 1.1], [0.6, 1.6], [1.2, 2.2]],
        bands_hz=[[1500, 3000], [3000, 6000], [6000, 12000]],
        flatness_db=rows,
    )


class CrashBalance:
    """Declared engineering ranking, never a substitute for listening approval."""

    def __init__(self, reference, rate):
        self.original = CrashObjective(reference, rate, 60)
        self.decay = BandDecayShapeLoss(reference, rate)
        shape = self.original.shape
        terminal = shape.target[:, -3:]
        plateau = (np.ptp(terminal, axis=1) < 1.5) & (
            shape.target.max(axis=1) - terminal.mean(axis=1) > 35
        )
        # Only quiet, flat terminal bands get a raised comparison floor.
        # Clearly decaying low rings are not treated as recording noise.
        self.tail_floor = np.where(plateau, terminal.mean(axis=1) + 10, -300)
        self.late = np.array([a >= 1.7 for a, b in shape.regions])
        self.specification = dict(
            spectral=self.original.specification,
            decay=self.decay.specification,
            weights=dict(
                mel=1, attack_mel=0.3, bloom=0.05, texture=0.15, shape_error_db=0.15
            ),
            acceptance="Inspect separate metrics, pitch, plots and held-out seeds; audition still required",
            version="crash-tail-floor-v3",
            terminal_plateau=dict(
                regions_seconds=[3, 6],
                max_range_db=1.5,
                below_peak_db=35,
                margin_db=10,
                apply_after_seconds=1.7,
                flagged_bands=np.flatnonzero(plateau).tolist(),
            ),
        )

    def components(self, audio):
        result = self.original.components(audio)
        result["bloom_unmasked"] = result["bloom"]
        shape = self.original.shape
        target, actual = shape.target.copy(), shape.db(shape.power(audio))
        target[:, self.late] = np.maximum(
            target[:, self.late], self.tail_floor[:, None]
        )
        actual[:, self.late] = np.maximum(
            actual[:, self.late], self.tail_floor[:, None]
        )
        error = (actual - target)[shape.active]
        rise = error[:, 2:9] - error[:, :1]
        result["bloom"] = float(np.sqrt(np.mean(error**2) + np.mean(rise**2)))
        result.update(self.decay.diagnostics(audio))
        result["score"] = (
            result["mel"]
            + 0.3 * result["attack_mel"]
            + 0.05 * result["bloom"]
            + 0.15 * result["texture"]
            + 0.15 * result["shape_error_db"]
        )
        return result


def plots(signals, rate, output):
    """Every plot uses fixed gains and identical analysis on source and synth."""
    bands = [
        (100, 300),
        (300, 700),
        (700, 1500),
        (1500, 3000),
        (3000, 6000),
        (6000, 16000),
    ]
    fig = make_subplots(
        rows=4,
        cols=2,
        subplot_titles=[
            "40–300 ms spectrum",
            "300–1000 ms spectrum",
            *[f"{lo}–{hi} Hz envelope" for lo, hi in bands],
        ],
    )
    for (name, audio), colour in zip(
        signals.items(), ["#eebc59", "#8d91a4", "#6bb9e9"]
    ):
        for col, (a, b) in enumerate(((0.04, 0.3), (0.3, 1)), 1):
            segment = audio[round(a * rate) : round(b * rate)]
            f, p = welch(segment, rate, nperseg=min(8192, len(segment)), nfft=16384)
            keep = (f >= 60) & (f <= 16000)
            fig.add_trace(
                go.Scatter(
                    x=f[keep],
                    y=10 * np.log10(np.maximum(p[keep], 1e-20)),
                    name=name,
                    legendgroup=name,
                    line_color=colour,
                    showlegend=col == 1,
                ),
                row=1,
                col=col,
            )
        f, t, z = stft(audio, rate, nperseg=4096, noverlap=3584)
        # A cropped render's zero-padded final window is not an audible burst.
        complete = t <= (len(audio) - 2048) / rate
        for i, (lo, hi) in enumerate(bands):
            power = np.sum(abs(z[(f >= lo) & (f < hi)]) ** 2, axis=0)
            fig.add_trace(
                go.Scatter(
                    x=t[complete],
                    y=10 * np.log10(np.maximum(power[complete], 1e-20)),
                    name=name,
                    legendgroup=name,
                    line_color=colour,
                    showlegend=False,
                ),
                row=2 + i // 2,
                col=1 + i % 2,
            )
    fig.update_xaxes(type="log", range=[np.log10(60), np.log10(16000)], row=1)
    fig.update_yaxes(range=[-100, -20], row=1)
    for row in range(2, 5):
        fig.update_xaxes(range=[0, 6], row=row)
        fig.update_yaxes(range=[-95, -10], row=row)
    fig.update_layout(
        template="plotly_dark",
        width=1450,
        height=1100,
        title="Crash: fixed-level spectrum and decay; no gain matching",
    )
    (output / "inspection.plotly.json").write_text(fig.to_json(), encoding="utf8")


def spectrograms(reference, audio, rate, output):
    """Reference-anchored colour and signed excess/missing energy, zero black."""
    spectra = [stft(x, rate, nperseg=4096, noverlap=3584) for x in (reference, audio)]
    f, t, _ = spectra[0]
    keep = (f >= 60) & (f < 16000)
    db = [20 * np.log10(np.maximum(abs(s[2][keep]), 1e-10)) for s in spectra]
    ceiling = float(db[0].max())
    floor = ceiling - 70
    fig = make_subplots(
        rows=3,
        cols=1,
        subplot_titles=[
            "Reference",
            "Candidate",
            "Difference: amber excess / blue missing",
        ],
    )
    for i, z in enumerate(db):
        fig.add_trace(
            go.Heatmap(
                x=t,
                y=np.log10(f[keep]),
                z=z,
                zmin=floor,
                zmax=ceiling,
                colorscale="Magma",
                showscale=False,
            ),
            row=i + 1,
            col=1,
        )
    diff = np.maximum(db[1], floor) - np.maximum(db[0], floor)
    fig.add_trace(
        go.Heatmap(
            x=t,
            y=np.log10(f[keep]),
            z=diff,
            zmin=-18,
            zmax=18,
            colorscale=[[0, "#218ad1"], [0.5, "#000000"], [1, "#efa740"]],
            showscale=False,
        ),
        row=3,
        col=1,
    )
    ticks = [100, 300, 1000, 3000, 10000]
    fig.update_yaxes(tickvals=np.log10(ticks), ticktext=[str(f) for f in ticks])
    fig.update_xaxes(range=[0, 4])
    fig.update_layout(template="plotly_dark", width=1450, height=1100)
    (output / "spectra.plotly.json").write_text(fig.to_json(), encoding="utf8")
