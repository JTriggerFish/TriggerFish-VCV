"""Show the measured bass sustain / noisy-attack tradeoff without auto gain."""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def balance_plot(loss, candidate):
    selected_bands = (0, 1, 5, 7)
    titles = [f"{loss.bands[i][0]:g}–{loss.bands[i][1]:g} Hz" for i in selected_bands]
    figure = make_subplots(
        rows=2, cols=2, subplot_titles=titles, shared_xaxes=True, shared_yaxes=True
    )
    candidate_power = loss.envelopes(candidate)
    for cell, band in enumerate(selected_bands):
        for name, powers, colour, dash in (
            ("Reference", loss.target_power, "#fcbf49", "solid"),
            ("Candidate", candidate_power, "#06d6a0", "dash"),
        ):
            figure.add_trace(
                go.Scatter(
                    x=loss.times,
                    y=10 * np.log10(np.maximum(powers[band], loss.floor)),
                    name=name,
                    showlegend=cell == 0,
                    line=dict(color=colour, dash=dash),
                ),
                row=cell // 2 + 1,
                col=cell % 2 + 1,
            )
    figure.update_xaxes(title="Seconds from onset", range=[0, 0.3])
    figure.update_yaxes(title="Filtered mean square, dBFS", range=[-90, -10])
    figure.update_layout(
        height=570, template="plotly_dark", margin=dict(l=70, r=30, t=45, b=40)
    )
    return (
        "<h2>Bass sustain and attack duration</h2>"
        "<p>Same causal band filters and power smoothing on both signals; "
        "fixed reference-derived floor, no level matching. Filter delay is "
        "measurement delay, not an inferred acoustic onset. These plots "
        "complement the spectral errors; they are not listening approval.</p>"
        + figure.to_html(full_html=False, include_plotlyjs=False)
    )
