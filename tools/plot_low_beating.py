"""Developer inspection: relative low-band pulse and modulation power, no audio gain edits."""

import json
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def main():
    directory = Path("build/crash-low-beating")
    audit = json.loads((directory / "audit.json").read_text())
    rows = {v["name"]: v["bands"][0] for v in audit["variants"]}
    figure = make_subplots(
        rows=2,
        cols=1,
        subplot_titles=[
            "90–180 Hz: pulsation after removing average exponential decay",
            "Modulation spectrum: slow beating versus rapid flutter",
        ],
    )
    for name, row, colour in [
        ("Reference", audit["reference"][0], "#f2bc53"),
        ("Current crash", rows["current"], "#64b6ec"),
        ("Close low pair experiment", rows["measured-pair-8"], "#c790d9"),
    ]:
        envelope = np.array(row["relative_envelope"])
        figure.add_trace(
            go.Scatter(
                x=0.5 + np.arange(len(envelope)) / 200,
                y=envelope,
                name=name,
                line_color=colour,
            ),
            row=1,
            col=1,
        )
        spectrum = np.array(row["spectrum"])
        figure.add_trace(
            go.Scatter(
                x=np.linspace(0, 100, len(spectrum)),
                y=10 * np.log10(np.maximum(spectrum, 1e-10)),
                name=name,
                showlegend=False,
                line_color=colour,
            ),
            row=2,
            col=1,
        )
    figure.update_xaxes(title_text="Time (s)", row=1, col=1)
    figure.update_xaxes(
        title_text="Pulses per second (Hz)", range=[0.5, 20], row=2, col=1
    )
    figure.update_yaxes(title_text="Relative envelope", row=1, col=1)
    figure.update_yaxes(
        title_text="Modulation power (dB/Hz)", range=[-70, 0], row=2, col=1
    )
    figure.update_layout(
        template="plotly_dark",
        height=1000,
        width=1450,
        title="Low-ring diagnostic — reference and synth audio levels are unchanged",
    )
    (directory / "beating.plotly.json").write_text(figure.to_json())


if __name__ == "__main__":
    main()
