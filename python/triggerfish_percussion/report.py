"""Self-contained Plotly A/B report for percussion references and renders."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import re

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .alignment import measure_alignment, shift_with_zeros
from .audio_io import AudioBuffer, resample_audio, write_wav
from .transforms import StftConfig, stft


@dataclass(frozen=True)
class ReportPair:
    label: str
    reference: AudioBuffer
    synthesis: AudioBuffer


@dataclass(frozen=True)
class AuditionClip:
    label: str
    audio: AudioBuffer


def write_comparison_report(
    pairs: tuple[ReportPair, ...],
    output: str | Path,
    title: str,
    auditions: tuple[AuditionClip, ...] = (),
) -> Path:
    if not pairs:
        raise ValueError("a comparison report requires at least one pair")
    path = Path(output)
    assets = path.with_name(path.stem + "-assets")
    assets.mkdir(parents=True, exist_ok=True)
    for pattern in ("*-reference.wav", "*-triggerfish.wav", "sweep-*.wav"):
        for stale in assets.glob(pattern):
            stale.unlink()
    sections: list[str] = []
    include_plotly = True
    for index, pair in enumerate(pairs):
        reference, synthesis = _aligned_pair(pair)
        audition_peak = max(
            float(np.max(np.abs(reference.samples))),
            float(np.max(np.abs(synthesis.samples))),
            1.0e-9,
        )
        common_gain = min(1.0, 0.9 / audition_peak)
        slug = _slug(f"{index + 1}-{pair.label}")
        reference_path = assets / f"{slug}-reference.wav"
        synthesis_path = assets / f"{slug}-triggerfish.wav"
        write_wav(reference_path, AudioBuffer(common_gain * reference.samples, 48000))
        write_wav(synthesis_path, AudioBuffer(common_gain * synthesis.samples, 48000))
        figure = _comparison_figure(reference, synthesis, pair.label)
        plot = figure.to_html(
            full_html=False,
            include_plotlyjs=True if include_plotly else False,
            config={"responsive": True, "displaylogo": False},
        )
        include_plotly = False
        relative_assets = assets.name
        sections.append(
            f"<section><h2>{escape(pair.label)}</h2>"
            f"<p>Both players use the same {common_gain:.4g} audition gain; "
            "analysis uses the untouched signals.</p>"
            "<div class='players'>"
            f"<label>Reference<audio controls preload='metadata' src='{relative_assets}/{reference_path.name}'></audio></label>"
            f"<label>TriggerFish<audio controls preload='metadata' src='{relative_assets}/{synthesis_path.name}'></audio></label>"
            f"</div>{plot}</section>"
        )
    audition_section = _write_auditions(auditions, assets)
    document = _html_document(title, audition_section + "".join(sections))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")
    return path


def _write_auditions(auditions: tuple[AuditionClip, ...], assets: Path) -> str:
    if not auditions:
        return ""
    players = []
    for index, clip in enumerate(auditions, 1):
        audio = resample_audio(clip.audio.mono(), 48000)
        peak = max(float(np.max(np.abs(audio.samples))), 1.0e-9)
        gain = min(1.0, 0.9 / peak)
        path = assets / f"sweep-{index}-{_slug(clip.label)}.wav"
        write_wav(path, AudioBuffer(gain * audio.samples, 48000))
        players.append(
            f"<label>{escape(clip.label)}"
            f"<audio controls preload='metadata' src='{assets.name}/{path.name}'></audio>"
            f"</label>"
        )
    return (
        "<section><h2>Current-model control sweeps</h2>"
        "<p>Even one-second hits; each file changes only its named control.</p>"
        f"<div class='players'>{''.join(players)}</div></section>"
    )


def _aligned_pair(pair: ReportPair) -> tuple[AudioBuffer, AudioBuffer]:
    reference = resample_audio(pair.reference.mono(), 48000)
    synthesis = resample_audio(pair.synthesis.mono(), 48000)
    length = min(reference.samples.size, synthesis.samples.size)
    first = reference.samples[:length]
    second = synthesis.samples[:length]
    alignment = measure_alignment(first, second, 48000)
    shifted = shift_with_zeros(second, alignment.candidate_lag_samples)
    return AudioBuffer(first, 48000), AudioBuffer(shifted, 48000)


def _comparison_figure(
    reference: AudioBuffer, synthesis: AudioBuffer, label: str
) -> go.Figure:
    seconds = np.arange(reference.samples.size) / reference.sample_rate
    stride = max(1, reference.samples.size // 12000)
    config = StftConfig(2048, 256, 4096)
    reference_stft = stft(reference.samples, reference.sample_rate, config)
    synthesis_stft = stft(synthesis.samples, synthesis.sample_rate, config)
    reference_db, synthesis_db = _common_magnitude_db(
        reference_stft.magnitude, synthesis_stft.magnitude, -100.0
    )
    figure = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=False,
        subplot_titles=(
            "Aligned waveform",
            "Reference spectrogram",
            "TriggerFish spectrogram",
        ),
        vertical_spacing=0.08,
    )
    figure.add_trace(
        go.Scattergl(
            x=seconds[::stride], y=reference.samples[::stride], name="Reference"
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scattergl(
            x=seconds[::stride], y=synthesis.samples[::stride], name="TriggerFish"
        ),
        row=1,
        col=1,
    )
    for row, result, magnitude_db in (
        (2, reference_stft, reference_db),
        (3, synthesis_stft, synthesis_db),
    ):
        selected = result.frequencies_hz >= 80.0
        frequencies = result.frequencies_hz[selected]
        times = result.times_seconds
        magnitude_db = magnitude_db[selected]
        frequency_stride = max(1, int(np.ceil(frequencies.size / 256)))
        time_stride = max(1, int(np.ceil(times.size / 400)))
        figure.add_trace(
            go.Heatmap(
                x=times[::time_stride],
                y=frequencies[::frequency_stride],
                z=magnitude_db[::frequency_stride, ::time_stride],
                zmin=-100,
                zmax=0,
                colorscale="Turbo",
                colorbar={"title": "dB"},
            ),
            row=row,
            col=1,
        )
        figure.update_yaxes(
            type="log", range=[np.log10(80), np.log10(22000)], row=row, col=1
        )
    figure.update_layout(title=label, height=1050, hovermode="x unified")
    figure.update_xaxes(title_text="Seconds", row=1, col=1)
    figure.update_xaxes(title_text="Seconds", row=3, col=1)
    return figure


def _common_magnitude_db(
    first: np.ndarray, second: np.ndarray, floor_db: float
) -> tuple[np.ndarray, np.ndarray]:
    peak = max(
        float(np.max(np.abs(first))),
        float(np.max(np.abs(second))),
        np.finfo(float).tiny,
    )
    floor = peak * 10.0 ** (floor_db / 20.0)

    def convert(values: np.ndarray) -> np.ndarray:
        return 20.0 * np.log10(np.maximum(np.abs(values), floor) / peak)

    return convert(first), convert(second)


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def _html_document(title: str, sections: str) -> str:
    return f"""<!doctype html><html><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width,initial-scale=1'>
<title>{escape(title)}</title><style>
body{{font-family:system-ui,sans-serif;max-width:1500px;margin:auto;padding:2rem;background:#111;color:#eee}}
section{{margin:2rem 0;padding:1rem;background:#181818;border-radius:.6rem}}
.players{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:1rem}}
label{{display:grid;gap:.4rem}}audio{{width:100%}}a{{color:#8cf}}
</style></head><body><h1>{escape(title)}</h1>
<p>Every comparison is real reference versus the current C++ model. There are no old-model renders.</p>
{sections}</body></html>"""
