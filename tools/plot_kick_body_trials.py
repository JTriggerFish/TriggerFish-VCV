"""Offline reference-aligned envelope plots for the body-hold experiment."""

import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.drum_balance_loss import DrumBalanceLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def main():
    root = Path(__file__).resolve().parents[1]
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset:]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        loss = DrumBalanceLoss(reference, rate)
        curves = [
            ("Reference", loss.target_power),
            ("Current", loss.envelopes(renderer.render(renderer.initial, 1.2))),
        ]
        for directory in (
            "kick-body-refinement",
            "kick-body-hold-20ms",
            "kick-body-hold-40ms",
        ):
            path = root / "build" / directory / "comparison.json"
            if path.exists():
                values = json.loads(path.read_text())["parameters"]
                curves.append((directory, loss.envelopes(renderer.render(values, 1.2))))
        figure = make_subplots(
            rows=2, cols=3, subplot_titles=[str(b) + " Hz" for b in loss.bands[:6]]
        )
        for index in range(6):
            for series, (name, power) in enumerate(curves):
                figure.add_trace(
                    go.Scatter(
                        x=loss.times,
                        y=10 * np.log10(np.maximum(power[index], loss.floor)),
                        name=name,
                        legendgroup=name,
                        showlegend=index == 0,
                        line=dict(
                            color=["orange", "grey", "cyan", "lime", "magenta"][series]
                        ),
                    ),
                    row=index // 3 + 1,
                    col=index % 3 + 1,
                )
        figure.update_xaxes(range=[0, 0.4], title="Seconds")
        figure.update_yaxes(range=[-80, -15], title="dBFS")
        figure.update_layout(template="plotly_dark", height=800, width=1400)
        figure.write_html(root / "build/kick-body-curves.html")
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
