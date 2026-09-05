"""Local-only reference/candidate reports; fixed amplitude and colour scales."""

import json
from html import escape
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.trajectory_fit_loss import TrajectoryLoss, features, to_db
from triggerfish_percussion.transforms import StftConfig, stft
from triggerfish_percussion.short_drum_fit_loss import ShortDrumLoss
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_publication import publish_html


def parameter_changes(directory):
    saved = json.loads((directory / "search.json").read_text(encoding="utf8"))
    rows = []
    for descriptor, initial in zip(
        saved["metadata"]["descriptors"], saved["metadata"]["values"]
    ):
        current = saved["parameters"][descriptor["key"]]
        if abs(current - initial) > 1e-6:
            rows.append(
                f'<tr><td>{descriptor["name"]}</td><td>{initial:.4g}</td>'
                f'<td>{current:.4g} {descriptor["unit"]}</td></tr>'
            )
    return (
        "<details><summary>Fitted controls — factory start → candidate</summary>"
        "<table><tr><th>Control</th><th>Start</th><th>Candidate</th></tr>"
        + "".join(rows)
        + "</table></details>"
    )


def write_report(directory: Path, kind="gong", assets="..", *, renderer):
    verify_candidate(renderer, directory)
    reference = read_wav(directory / "reference.wav").mono()
    candidate = read_wav(directory / "candidate.wav").mono()
    baseline = read_wav(directory / "baseline.wav").mono()
    saved = json.loads((directory / "search.json").read_text(encoding="utf8"))
    source = saved["metadata"]["reference"]
    cell = source.get("cell", {})
    identity = (
        f'<p>Target: {escape(source["name"])} · '
        f'{escape(str(cell.get("articulation", "")))} · '
        f'velocity {cell.get("velocity", "")} · variation {cell.get("repeat", "")}. '
        f'Source gain {source.get("referenceGainDb", 0):g} dB; '
        f'onset {1000 * cell.get("onset_seconds", 0):.3f} ms.</p>'
    )
    loss_type = ShortDrumLoss if kind == "kick" else TrajectoryLoss
    loss = loss_type(reference.samples, reference.sample_rate)
    metrics = {
        name: loss.diagnostics(audio.samples)
        for name, audio in [("baseline", baseline), ("candidate", candidate)]
    }
    (directory / "comparison.json").write_text(json.dumps(metrics, indent=2))
    transforms = [
        stft(
            audio.samples,
            audio.sample_rate,
            StftConfig(
                8192 if kind == "kick" else 4096, 256 if kind == "kick" else 512
            ),
        )
        for audio in (reference, candidate)
    ]
    maximum = 10 * np.log10(max(np.max(transforms[0].power), 1e-30))
    figure = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        subplot_titles=["Reference", "Candidate — exact workbench renderer"],
    )
    for row, value in enumerate(transforms, 1):
        selected = (value.frequencies_hz >= (20 if kind == "kick" else 40)) & (
            value.frequencies_hz <= 16000
        )
        figure.add_trace(
            go.Heatmap(
                x=value.times_seconds,
                y=value.frequencies_hz[selected],
                z=np.clip(
                    to_db(value.power[selected], 1e-20), maximum - 80, maximum
                ).tolist(),
                zmin=maximum - 80,
                zmax=maximum,
                coloraxis="coloraxis",
            ),
            row=row,
            col=1,
        )
        figure.update_yaxes(type="log", title="Hz", row=row, col=1)
    figure.update_xaxes(title="Seconds from onset", row=2, col=1)
    figure.update_layout(
        height=650,
        template="plotly_dark",
        coloraxis=dict(
            colorscale="Magma",
            cmin=maximum - 80,
            cmax=maximum,
            colorbar=dict(title="dBFS"),
        ),
        margin=dict(l=65, r=60, t=45, b=40),
    )
    trajectory = make_subplots(rows=1, cols=1)
    colours = ["#fcbf49", "#ef476f", "#06d6a0", "#118ab2", "#bc89ff"]
    bands = [(40, 300), (300, 1000), (1000, 3000), (3000, 8000), (8000, 16000)]
    if kind == "kick":
        bands = [(20, 120), (120, 300), (300, 1000), (1000, 3000), (3000, 16000)]
    for name, audio, dash in [
        ("Reference", reference, "solid"),
        ("Candidate", candidate, "dash"),
    ]:
        value = features(audio.samples, audio.sample_rate) if kind != "kick" else None
        if kind == "kick":
            value = stft(audio.samples, audio.sample_rate, StftConfig(2048, 128))
            centres, powers, times = (
                value.frequencies_hz,
                value.power,
                value.times_seconds,
            )
        else:
            centres, powers, times = value.band_centres, value.bands, value.times
        for (low, high), colour in zip(bands, colours):
            selected = (centres >= low) & (centres < high)
            trajectory.add_trace(
                go.Scatter(
                    x=times,
                    y=to_db(powers[selected].sum(axis=0), 1e-14),
                    name=f"{name} {low}–{high} Hz",
                    line=dict(color=colour, dash=dash),
                )
            )
    trajectory.update_layout(
        height=400,
        template="plotly_dark",
        xaxis_title="Seconds",
        yaxis_title="Band power, dBFS",
        yaxis_range=[-100, -10],
    )
    keys = (
        (
            "attack_spectrum_rmse_db",
            "body_spectrum_rmse_db",
            "low_spectrum_rmse_db",
            "envelope_rmse_db",
        )
        if kind == "kick"
        else ("band_rmse_db", "ridge_rmse_db")
    )
    headings = (
        ("Attack spectrum", "Body spectrum", "Low spectrum", "RMS envelope")
        if kind == "kick"
        else ("Band power", "Low spectrum")
    )
    table = "<p>Candidate errors against reference (dB; not a perceptual score)</p><table><tr><th>Region (s)</th>"
    table += "".join(f"<th>{name}</th>" for name in headings) + "</tr>"
    for after in metrics["candidate"]["regions"]:
        table += f'<tr><td>{after["seconds"]}</td>'
        table += "".join(f"<td>{after[key]:.1f}</td>" for key in keys) + "</tr>"
    table += "</table>"
    html = (
        Path(__file__)
        .with_name("workbench_fit_report_shell.html")
        .read_text(encoding="utf8")
        .replace("../lookahead_", f"{assets}/lookahead_")
        + identity
        + table
        + parameter_changes(directory)
        + figure.to_html(
            full_html=False, include_plotlyjs=f"{assets}/vendor/plotly.min.js"
        )
        + trajectory.to_html(full_html=False, include_plotlyjs=False)
    )
    if kind == "kick":
        from kick_fit_plots import detail_plot

        html += detail_plot(reference, candidate)
        sequence_buttons = []
        for name, label in (
            ("hits-0.5s-strength-0.5", "Quarter notes · medium"),
            ("hits-0.125s-strength-0.5", "Fast retriggers · medium"),
            ("hits-0.125s-strength-1", "Fast retriggers · full strength"),
        ):
            if (directory / f"{name}.wav").exists():
                sequence_buttons.append(f'<button data-audio="{name}">{label}</button>')
        if sequence_buttons:
            html += "<h2>Playing checks · same candidate</h2><p>"
            html += " ".join(sequence_buttons) + "</p>"
        html = html.replace("Gong", "Acoustic kick").replace(
            "One mallet strike", "One centre beater strike"
        )
    publish_html(directory, html)


if __name__ == "__main__":
    import sys

    raise SystemExit(
        "Use dev.ps1 fit-kick-start / fit-gong-start: reports require renderer verification"
    )
