"""Ablate turbulence mechanisms and measure tail-envelope modulation.

The dB trend is removed for analysis only; WAVs retain their actual levels.
Different exchange settings consume different random streams, so compare several
seeds statistically, not waveform differences between nominally identical seeds.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, sosfilt, periodogram
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.utils import PlotlyJSONEncoder

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def envelope(samples, rate, band):
    filtered = sosfilt(
        butter(2, band, btype="bandpass", fs=rate, output="sos"), samples
    )
    # Decimate power before smoothing to keep this diagnostic inexpensive.
    hop = round(rate * 0.004)
    frames = len(samples) // hop
    power = (filtered[: frames * hop] ** 2).reshape(frames, hop).mean(axis=1)
    power = gaussian_filter1d(power, 3)
    return np.arange(frames) * hop / rate, 10 * np.log10(np.maximum(power, 1e-20))


def measure(times, db, mask):
    trend = gaussian_filter1d(db, 0.3 / (times[1] - times[0]))
    fluctuation = (db - trend)[mask]
    frequencies, power = periodogram(fluctuation, fs=1 / (times[1] - times[0]))
    selected = (frequencies >= 1) & (frequencies <= 30)
    peak = np.flatnonzero(selected)[np.argmax(power[selected])]
    return dict(
        ripple_rms_db=float(np.std(fluctuation)),
        strongest_modulation_hz=float(frequencies[peak]),
        strongest_bin_fraction=float(power[peak] / max(1e-30, power[selected].sum())),
    )


def run(args):
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        parameters = dict(renderer.initial)
        if args.cascade:
            parameters["bloom_rate"] = args.cascade
        args.output.mkdir(parents=True, exist_ok=True)
        seconds = 10 if args.target == "crash" else 8
        reference = aligned_reference(renderer, seconds)
        rate = renderer.sample_rate
        bands = [(100, 700), (700, 3000), (3000, 16000)]
        reference_curves = [envelope(reference, rate, band) for band in bands]
        masks = [
            (t >= 1) & (t <= 5) & (db > db.max() - 45) for t, db in reference_curves
        ]
        figure = make_subplots(
            rows=3, cols=1, subplot_titles=[f"{a}–{b} Hz" for a, b in bands]
        )
        report = dict(
            parameters=parameters, reference=renderer.metadata["reference"], variants={}
        )
        variants = {
            "reference": None,
            "relaxed": {},
            "classic": {"field_relaxed_turbulence": 0},
            "exchange-off": {"field_exchange": 0},
            "half-phase": {
                "field_phase_bandwidth": parameters["field_phase_bandwidth"] * 0.5
            },
            "phase-off": {"field_phase_bandwidth": 0},
            "full-density": {"field_satellite_density": 1},
        }
        seed = renderer.metadata["event"]["seed"]
        for label, changes in variants.items():
            rows = []
            for offset in ((0,) if changes is None else (0, 13001, 14009)):
                patch = (
                    dict(parameters, field_relaxed_turbulence=1, **(changes or {}))
                    if changes is None or "field_relaxed_turbulence" not in changes
                    else dict(parameters, **changes)
                )
                samples = (
                    reference
                    if changes is None
                    else renderer.render(patch, seconds, (seed + offset) & 0xFFFFFFFF)
                )
                curves = [envelope(samples, rate, band) for band in bands]
                rows.append(
                    dict(
                        seed=(seed + offset) & 0xFFFFFFFF,
                        bands=[
                            measure(t, db, mask) for (t, db), mask in zip(curves, masks)
                        ],
                    )
                )
                if offset == 0:
                    write_wav(args.output / f"{label}.wav", AudioBuffer(samples, rate))
                    for index, (t, db) in enumerate(curves):
                        figure.add_trace(
                            go.Scatter(
                                x=t,
                                y=db,
                                name=label,
                                legendgroup=label,
                                showlegend=index == 0,
                            ),
                            row=index + 1,
                            col=1,
                        )
            report["variants"][label] = rows
        figure.update_xaxes(title_text="Seconds", range=[0, 6])
        figure.update_yaxes(title_text="Band power (dB)", range=[-100, -10])
        figure.update_layout(
            template="plotly_dark",
            width=1450,
            height=1050,
            title="Tail modulation ablation — unchanged WAV levels",
        )
        (args.output / "modulation.plotly.json").write_text(
            json.dumps(figure.to_plotly_json(), cls=PlotlyJSONEncoder)
        )
        (args.output / "modulation.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report["variants"]), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cascade", type=float)
    run(parser.parse_args())
