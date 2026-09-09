"""Controlled exact-Wasm comparison; no fitting, gain matching or preset writes.

Hold contact, painted modes, turbulence, damping and observation fixed. Sweep
only the declared transfer mode/strength/nonlinearity. Save replayable workbench
fits and fixed-level WAVs; compare band evolution and auraloss over three seeds.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.utils import PlotlyJSONEncoder
import torch

from check_metal_tail_modulation import envelope
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer

BANDS = [(80, 700), (700, 3000), (3000, 16000)]
COLOURS = ["#ab63fa", "#ff6692", "#fecb52", "#00cc96", "#19d3f3", "#ff97ff"]


def measures(samples, rate, reference_curves):
    curves = [envelope(samples, rate, band) for band in BANDS]
    rows = []
    for (time, db), (_, reference_db) in zip(curves, reference_curves):
        active = reference_db > reference_db.max() - 45
        rows.append(
            dict(
                peak_time=float(time[np.argmax(db)]),
                level_error_db=float(
                    np.sqrt(np.mean((db[active] - reference_db[active]) ** 2))
                ),
                tail_relative_db=float(
                    np.mean(db[(time >= 1) & (time < 3)]) - db.max()
                ),
            )
        )
    return rows, curves


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        args.output.mkdir(parents=True, exist_ok=True)
        seconds = 10 if args.target == "crash" else 8
        reference = aligned_reference(renderer, seconds)
        rate = renderer.sample_rate
        reference_curves = [envelope(reference, rate, band) for band in BANDS]
        perceptual = AuralossMel(reference, rate)
        base = dict(renderer.initial)
        variants = {
            "current-diffusion": {},
        }
        for strength in (0.25, 1, 4, 16):
            variants[f"quadratic-{strength:g}"] = dict(
                bloom_rate=strength,
                bloom_energy_acceleration=1,
            )
        report = dict(
            metadata=renderer.metadata,
            baseline_parameters=base,
            variants={},
            audio_normalized=False,
            fitting_performed=False,
            listening_approved=False,
        )
        figure = make_subplots(
            rows=3, cols=1, subplot_titles=[f"{a}–{b} Hz" for a, b in BANDS]
        )
        write_wav(args.output / "reference.wav", AudioBuffer(reference, rate))
        for index, (time, db) in enumerate(reference_curves):
            figure.add_trace(
                go.Scatter(
                    x=time,
                    y=db,
                    name="Reference",
                    legendgroup="reference",
                    line=dict(color="white", width=3),
                    showlegend=index == 0,
                ),
                row=index + 1,
                col=1,
            )
        seed = renderer.metadata["event"]["seed"]
        for variant_index, (label, overrides) in enumerate(variants.items()):
            parameters = dict(base, **overrides)
            rows = []
            for offset in (0, 17011, 18013):
                samples = renderer.render(
                    parameters, seconds, (seed + offset) & 0xFFFFFFFF
                )
                if not np.isfinite(samples).all():
                    raise ValueError(f"Nonfinite {label} render")
                bands, curves = measures(samples, rate, reference_curves)
                rows.append(
                    dict(
                        seed=(seed + offset) & 0xFFFFFFFF,
                        bands=bands,
                        mel=perceptual.score(samples),
                        peak_db=float(20 * np.log10(max(1e-15, abs(samples).max()))),
                    )
                )
                if not offset:
                    write_wav(args.output / f"{label}.wav", AudioBuffer(samples, rate))
                    snapshot = renderer.request(
                        command="snapshot", parameters=parameters, name=label
                    )["fit"]
                    replay = renderer.decode(
                        renderer.request(
                            command="renderSnapshot", fit=snapshot, seconds=seconds
                        )["pcm"]
                    )
                    if not np.array_equal(samples, replay):
                        raise ValueError(f"Snapshot does not exactly replay {label}")
                    (args.output / f"{label}.fit.json").write_text(
                        json.dumps(snapshot, indent=2)
                    )
                    for index, (time, db) in enumerate(curves):
                        figure.add_trace(
                            go.Scatter(
                                x=time,
                                y=db,
                                name=label,
                                legendgroup=label,
                                line=dict(color=COLOURS[variant_index]),
                                showlegend=index == 0,
                            ),
                            row=index + 1,
                            col=1,
                        )
            repeated = renderer.decode(
                renderer.request(
                    command="renderSequence",
                    parameters=parameters,
                    seconds=8,
                    hits=[dict(time=i * 0.5, seed=seed + i) for i in range(8)],
                )["pcm"]
            )
            if not np.isfinite(repeated).all():
                raise ValueError(f"Nonfinite repeated-hit {label} render")
            write_wav(
                args.output / f"{label}-quarters.wav", AudioBuffer(repeated, rate)
            )
            report["variants"][label] = dict(
                overrides=overrides,
                seeds=rows,
                exact_snapshot_replay=True,
                repeated_peak_db=float(20 * np.log10(max(1e-15, abs(repeated).max()))),
            )
            print(label, json.dumps(rows), flush=True)
        figure.update_xaxes(title_text="Seconds", range=[0, 6])
        figure.update_yaxes(title_text="Band power (dB)", range=[-100, -10])
        figure.update_layout(
            template="plotly_dark",
            width=1450,
            height=1050,
            title=f"{args.target}: controlled transport comparison — unchanged levels",
        )
        (args.output / "comparison.plotly.json").write_text(
            json.dumps(figure.to_plotly_json(), cls=PlotlyJSONEncoder)
        )
        (args.output / "comparison.json").write_text(json.dumps(report, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
