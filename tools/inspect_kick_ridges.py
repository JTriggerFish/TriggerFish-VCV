"""Inspect actual narrow spectral excesses, not only broad-band totals."""

import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import find_peaks

from triggerfish_percussion.transforms import StftConfig, stft
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "build/kick-ridges"
    out.mkdir(exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        ref = renderer.reference[onset:]
        ref = np.pad(ref, (0, round(0.8 * rate) - len(ref)))
        values = renderer.initial
        if os.environ.get("TF_KICK_CANDIDATE"):
            fit = json.loads(
                Path(os.environ["TF_KICK_CANDIDATE"]).read_text(encoding="utf8")
            )
            values = {
                k: v
                for n in fit["instrument"]["nodes"]
                for k, v in n["parameters"].items()
            }
        synth = renderer.render(values, 0.8)
        signals = [("Reference", ref), ("Kick", synth)]
        for key in ("thump_level", "resonance_level"):
            p = dict(values, thump_level=0, resonance_level=0, contact_level=0)
            p[key] = values[key]
            signals.append((key, renderer.render(p, 0.8)))
        transforms = [stft(x, rate, StftConfig(8192, 128)) for _, x in signals]
        freq, times = transforms[0].frequencies_hz, transforms[0].times_seconds
        fig = make_subplots(
            rows=3,
            cols=2,
            subplot_titles=[
                "Reference",
                "Kick",
                "0–80 ms",
                "80–160 ms",
                "160–260 ms",
                "260–400 ms",
            ],
        )
        ceiling = 10 * np.log10(transforms[0].power.max())
        for i in range(2):
            t = transforms[i]
            bins = (freq >= 20) & (freq <= 2000)
            fig.add_trace(
                go.Heatmap(
                    x=times,
                    y=freq[bins],
                    z=10 * np.log10(np.maximum(t.power[bins], 1e-20)),
                    zmin=ceiling - 70,
                    zmax=ceiling,
                    showscale=False,
                    colorscale="Magma",
                ),
                row=1,
                col=i + 1,
            )
            fig.update_yaxes(
                type="log", range=[np.log10(20), np.log10(2000)], row=1, col=i + 1
            )
            fig.update_xaxes(range=[0, 0.4], row=1, col=i + 1)
        rows = []
        for index, (start, end) in enumerate(
            ((0, 0.08), (0.08, 0.16), (0.16, 0.26), (0.26, 0.4))
        ):
            power = [
                t.power[:, (times >= start) & (times < end)].mean(axis=1)
                for t in transforms
            ]
            db = [10 * np.log10(np.maximum(p, 1e-20)) for p in power]
            peaks, _ = find_peaks(db[1], prominence=3)
            peaks = [
                i for i in peaks if 20 < freq[i] < 2000 and db[1][i] > ceiling - 60
            ]
            peaks = sorted(peaks, key=lambda i: db[1][i] - db[0][i], reverse=True)
            rows.append(
                dict(
                    region=[start, end],
                    excess_peaks=[
                        dict(
                            hz=float(freq[i]),
                            excess_db=float(db[1][i] - db[0][i]),
                            reference_db=float(db[0][i]),
                            synth_db=float(db[1][i]),
                        )
                        for i in peaks[:8]
                    ],
                )
            )
            reference_peaks, _ = find_peaks(db[0], prominence=2)
            rows[-1]["reference_bass_peaks"] = [
                dict(hz=float(freq[i]), db=float(db[0][i]))
                for i in reference_peaks
                if 20 < freq[i] < 350 and db[0][i] > ceiling - 65
            ]
            for j, (name, _) in enumerate(signals):
                fig.add_trace(
                    go.Scatter(
                        x=freq,
                        y=db[j],
                        name=name,
                        legendgroup=name,
                        showlegend=index == 0,
                        line=dict(color=["orange", "cyan", "pink", "lime"][j]),
                    ),
                    row=2 + index // 2,
                    col=1 + index % 2,
                )
            fig.update_xaxes(
                type="log",
                range=[np.log10(20), np.log10(2000)],
                row=2 + index // 2,
                col=1 + index % 2,
            )
            fig.update_yaxes(
                range=[ceiling - 80, ceiling + 5], row=2 + index // 2, col=1 + index % 2
            )
        fig.update_layout(template="plotly_dark", width=1450, height=1000)
        fig.write_html(out / "ridges.html")
        report = dict(parameters=values, metadata=renderer.metadata, ridges=rows)
        (out / "ridges.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(rows, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
