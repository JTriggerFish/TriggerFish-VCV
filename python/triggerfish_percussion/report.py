"""Self-contained Plotly A/B report for percussion references and renders."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import re

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .alignment import detect_impact_onset, measure_alignment, shift_with_zeros
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
    acceptance: dict[str, object] | None = None,
    causal_fit: dict[str, object] | None = None,
    audition_gain: float = 1.0,
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
    common_gain = _safe_report_gain(audition_gain, pairs, auditions)
    gain_db = 20.0 * np.log10(common_gain)
    cache_token = f"gain-{round(common_gain * 1.0e6)}"
    for index, pair in enumerate(pairs):
        reference, synthesis = _aligned_pair(pair)
        slug = _slug(f"{index + 1}-{pair.label}")
        reference_path = assets / f"{slug}-reference.wav"
        synthesis_path = assets / f"{slug}-triggerfish.wav"
        write_wav(reference_path, AudioBuffer(common_gain * reference.samples, 48000))
        write_wav(synthesis_path, AudioBuffer(common_gain * synthesis.samples, 48000))
        initial_reference, initial_synthesis = _onset_pair(reference, synthesis, 0.100)
        initial_reference_path = assets / f"{slug}-first-100ms-reference.wav"
        initial_synthesis_path = assets / f"{slug}-first-100ms-triggerfish.wav"
        write_wav(
            initial_reference_path,
            AudioBuffer(common_gain * initial_reference.samples, 48000),
        )
        write_wav(
            initial_synthesis_path,
            AudioBuffer(common_gain * initial_synthesis.samples, 48000),
        )
        repeated_reference_path = assets / f"{slug}-first-100ms-repeat-reference.wav"
        repeated_synthesis_path = assets / f"{slug}-first-100ms-repeat-triggerfish.wav"
        write_wav(
            repeated_reference_path,
            _repeat_with_gap(initial_reference, common_gain),
        )
        write_wav(
            repeated_synthesis_path,
            _repeat_with_gap(initial_synthesis, common_gain),
        )
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
            f"<p>Every report clip uses one source-grid audition gain "
            f"({gain_db:+.1f} dB); analysis uses the "
            "untouched signals. Short clips begin at the detected impact.</p>"
            "<h3>First 100 ms only</h3>"
            "<div class='players'>"
            f"<label>Reference — 100 ms<audio controls preload='metadata' src='{relative_assets}/{initial_reference_path.name}?v={cache_token}'></audio></label>"
            f"<label>TriggerFish — 100 ms<audio controls preload='metadata' src='{relative_assets}/{initial_synthesis_path.name}?v={cache_token}'></audio></label>"
            "</div><h3>First 100 ms repeated four times</h3>"
            "<div class='players'>"
            f"<label>Reference — repeated<audio controls preload='metadata' src='{relative_assets}/{repeated_reference_path.name}?v={cache_token}'></audio></label>"
            f"<label>TriggerFish — repeated<audio controls preload='metadata' src='{relative_assets}/{repeated_synthesis_path.name}?v={cache_token}'></audio></label>"
            "</div><h3>Full recording</h3><div class='players'>"
            f"<label>Reference<audio controls preload='metadata' src='{relative_assets}/{reference_path.name}?v={cache_token}'></audio></label>"
            f"<label>TriggerFish<audio controls preload='metadata' src='{relative_assets}/{synthesis_path.name}?v={cache_token}'></audio></label>"
            f"</div>{plot}</section>"
        )
    audition_section = _write_auditions(auditions, assets, common_gain, cache_token)
    gate_section = _acceptance_section(acceptance)
    causal_section = _causal_fit_section(causal_fit)
    document = _html_document(
        title, gate_section + causal_section + audition_section + "".join(sections)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")
    return path


def velocity_grid_audition_gain(
    references: tuple[tuple[float, AudioBuffer], ...],
    target_average_peak_dbfs: float = -12.0,
    headroom_dbfs: float = -1.0,
) -> float:
    """Choose one gain from equally weighted source-velocity peak levels."""
    if not references:
        raise ValueError("audition normalization needs source velocity references")
    levels: dict[float, list[float]] = {}
    maximum_peak = 0.0
    for velocity, audio in references:
        peak = float(np.max(np.abs(audio.mono().samples)))
        if peak <= 0.0 or not np.isfinite(peak):
            continue
        maximum_peak = max(maximum_peak, peak)
        levels.setdefault(float(velocity), []).append(20.0 * np.log10(peak))
    if not levels:
        raise ValueError("audition normalization references are silent")
    velocity_levels = [float(np.mean(values)) for values in levels.values()]
    requested = 10.0 ** (
        (target_average_peak_dbfs - float(np.mean(velocity_levels))) / 20.0
    )
    headroom = 10.0 ** (headroom_dbfs / 20.0) / maximum_peak
    return float(min(requested, headroom))


def _acceptance_section(diagnostics: dict[str, object] | None) -> str:
    if diagnostics is None:
        return ""
    passed = bool(diagnostics["passed"])
    rows = []
    texture_rows = []
    for cell in diagnostics["cells"]:
        for region in cell["regions"]:
            state = "pass" if region["passed"] else "FAIL"
            rows.append(
                "<tr>"
                f"<td>{escape(cell['cell'])}</td><td>{escape(region['region'])}</td>"
                f"<td>{region['level_delta_db']:+.2f}</td>"
                f"<td>{region['erb_shape_rmse_db']:.2f}</td>"
                f"<td>{region['maximum_absolute_erb_delta_db']:.2f}</td>"
                f"<td>{state}</td></tr>"
            )
        trajectory = cell.get("trajectory")
        if trajectory is not None:
            state = "pass" if trajectory["passed"] else "FAIL"
            rows.append(
                "<tr>"
                f"<td>{escape(cell['cell'])}</td><td>time trajectory</td>"
                "<td>—</td>"
                f"<td>{trajectory['rmse_db']:.2f}</td>"
                f"<td>{trajectory['maximum_absolute_delta_db']:.2f}</td>"
                f"<td>{state}</td></tr>"
            )
        stability = cell.get("stability")
        if stability is not None:
            state = "pass" if stability["passed"] else "FAIL"
            rows.append(
                "<tr>"
                f"<td>{escape(cell['cell'])}</td><td>late regrowth</td>"
                f"<td>{stability['maximum_late_regrowth_db']:+.2f}</td>"
                "<td>—</td><td>—</td>"
                f"<td>{state}</td></tr>"
            )
        texture = cell.get("initial_texture")
        if texture is not None:
            state = "pass" if texture["passed"] else "FAIL"
            texture_rows.append(
                "<tr>"
                f"<td>{escape(cell['cell'])}</td>"
                f"<td>{texture['fine_spectrum_rmse_db']:.2f}</td>"
                f"<td>{texture['centroid_rmse_octaves']:.3f}</td>"
                f"<td>{texture['rolloff_rmse_octaves']:.3f}</td>"
                f"<td>{texture['flatness_rmse_db']:.2f}</td>"
                f"<td>{texture['crest_rmse_db']:.2f}</td>"
                f"<td>{texture['ridge_ratio_absolute_error']:.3f}</td>"
                f"<td>{state}</td></tr>"
            )
    modes = diagnostics.get("persistent_modes")
    mode_summary = ""
    if modes is not None:
        state = "pass" if modes["passed"] else "FAIL"
        mode_summary = (
            "<p><strong>Persistent modes:</strong> "
            f"{modes['matched_count']}/{modes['reference_count']} reference modes "
            f"matched; frequency RMSE {modes['frequency_rmse_cents']:.1f} cents; "
            f"gate {state}.</p>"
        )
    verdict = "ACCEPTED" if passed else "REJECTED"
    texture_table = ""
    if texture_rows:
        texture_table = (
            "<h3>Initial 100 ms brightness and texture</h3>"
            "<table><thead><tr><th>Cell</th><th>Fine spectrum RMSE dB</th>"
            "<th>Centroid RMSE oct</th>"
            "<th>Rolloff RMSE oct</th><th>Flatness RMSE dB</th>"
            "<th>Crest RMSE dB</th><th>Ridge-ratio error</th><th>Gate</th>"
            f"</tr></thead><tbody>{''.join(texture_rows)}</tbody></table>"
        )
    return (
        f"<section class='verdict'><h2>Calibration gates: {verdict}</h2>"
        "<p>A failed region cannot be waived by the aggregate optimizer score.</p>"
        f"{mode_summary}"
        "<table><thead><tr><th>Cell</th><th>Region</th><th>Level Δ dB</th>"
        "<th>ERB shape RMSE dB</th><th>Worst ERB Δ dB</th><th>Gate</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table>{texture_table}</section>"
    )


def _causal_fit_section(diagnostics: dict[str, object] | None) -> str:
    if diagnostics is None:
        return ""
    rows = []
    for stage in diagnostics["stages"]:
        key = f"{stage['end_seconds']:.3f}"
        baseline = stage["baseline_prefix_losses"][key]
        candidate = stage["candidate_prefix_losses"][key]
        update = "updated" if stage["accepted"] else "kept baseline"
        quality = stage.get("candidate_absolute_quality", stage.get("absolute_quality"))
        if quality is None:
            quality_state = "not recorded"
            envelope = spectral = spectral_p95 = float("nan")
            texture_summary = "not recorded"
        else:
            cells = quality["cells"]
            envelope = max(cell["envelope_rmse_db"] for cell in cells)
            spectral = max(cell["spectral_rmse_db"] for cell in cells)
            spectral_p95 = max(cell["spectral_p95_absolute_db"] for cell in cells)
            texture_summary = " / ".join(
                (
                    f"C {max(cell.get('centroid_rmse_octaves', 0.0) for cell in cells):.2f} oct",
                    f"S {max(cell.get('fine_spectrum_rmse_db', 0.0) for cell in cells):.1f} dB",
                    f"R {max(cell.get('rolloff_rmse_octaves', 0.0) for cell in cells):.2f} oct",
                    f"F {max(cell.get('flatness_rmse_db', 0.0) for cell in cells):.1f} dB",
                    f"P {max(cell.get('crest_rmse_db', 0.0) for cell in cells):.1f} dB",
                    f"ridge {max(cell.get('ridge_ratio_absolute_error', 0.0) for cell in cells):.2f}",
                )
            )
            if not quality["required"]:
                quality_state = "initialization candidate"
            else:
                quality_state = "PASS" if quality["passed"] else "FAIL"
        cumulative = stage.get("candidate_cumulative_absolute_quality", ())
        cumulative_state = (
            "PASS"
            if cumulative and all(item["passed"] for item in cumulative)
            else ("FAIL" if cumulative else "—")
        )
        retry = " (joint retry)" if stage.get("retry") else ""
        rows.append(
            "<tr>"
            f"<td>{escape(stage['stage'] + retry)}</td>"
            f"<td>{stage['end_seconds']:g}</td>"
            f"<td>{baseline:.6g}</td><td>{candidate:.6g}</td>"
            f"<td>{stage['worst_protected_loss_ratio']:.4f}</td>"
            f"<td>{stage['worst_current_loss_ratio']:.4f}</td>"
            f"<td>{update}</td>"
            f"<td>{envelope:.2f}</td><td>{spectral:.2f}</td>"
            f"<td>{spectral_p95:.2f}</td><td>{texture_summary}</td>"
            f"<td>{quality_state}</td>"
            f"<td>{cumulative_state}</td></tr>"
        )
    complete = diagnostics.get("completed")
    blocked = diagnostics.get("blocked_stage")
    if complete is True:
        verdict = "All active-prefix gates passed."
    elif blocked:
        verdict = (
            f"Fitting stopped at {escape(str(blocked))}: absolute prefix gate failed."
        )
    else:
        verdict = "Completion state was not recorded."
    policy = diagnostics.get("policy", "strict")
    configuration = diagnostics.get("policy_configuration", {})
    default_tolerance = 0.35 if policy == "first-100ms-tradeoff" else 0.025
    tolerance = 100.0 * float(
        configuration.get("protected_acceptance_tolerance", default_tolerance)
    )
    if policy == "first-100ms-tradeoff":
        policy_description = (
            "The 100 ms descriptor owns the objective; nested 4 and 15 ms "
            "descriptors have reduced weights so their attack features are not "
            f"triple-counted. Earlier per-cell loss may rise by at most {tolerance:g}%. "
            "Absolute level, spectrum, brightness, and texture errors remain hard "
            "acceptance barriers."
        )
    else:
        policy_description = (
            "Earlier per-cell prefixes may not exceed their retained "
            f"{tolerance:g}% loss ceiling. Absolute errors are hard acceptance gates."
        )
    return (
        "<section><h2>Cumulative causal fit</h2>"
        f"<p><strong>{verdict}</strong></p>"
        "<p>Every row uses the complete render from sample zero. "
        f"Policy: <code>{escape(str(policy))}</code>. {policy_description} "
        "The tabulated absolute errors are worst-cell values.</p>"
        "<table><thead><tr><th>Stage</th><th>Prefix s</th>"
        "<th>Baseline loss</th><th>Candidate loss</th>"
        "<th>Worst prior ratio</th><th>Worst current ratio</th><th>Optimizer</th>"
        "<th>Candidate envelope RMSE dB</th><th>Candidate spectral RMSE dB</th>"
        "<th>Candidate spectral p95 |Δ| dB</th><th>Worst texture errors</th>"
        "<th>Active-prefix gate</th>"
        "<th>All-prefix gates</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table></section>"
    )


def _write_auditions(
    auditions: tuple[AuditionClip, ...],
    assets: Path,
    gain: float,
    cache_token: str,
) -> str:
    if not auditions:
        return ""
    players = []
    for index, clip in enumerate(auditions, 1):
        audio = resample_audio(clip.audio.mono(), 48000)
        path = assets / f"sweep-{index}-{_slug(clip.label)}.wav"
        write_wav(path, AudioBuffer(gain * audio.samples, 48000))
        players.append(
            f"<label>{escape(clip.label)}"
            f"<audio controls preload='metadata' src='{assets.name}/{path.name}?v={cache_token}'></audio>"
            f"</label>"
        )
    return (
        "<section><h2>Current-model control sweeps</h2>"
        "<p>Even one-second hits; each file changes only its named control.</p>"
        f"<div class='players'>{''.join(players)}</div></section>"
    )


def _safe_report_gain(
    requested: float,
    pairs: tuple[ReportPair, ...],
    auditions: tuple[AuditionClip, ...],
) -> float:
    if not np.isfinite(requested) or requested <= 0.0:
        raise ValueError("report audition gain must be finite and positive")
    peak = max(
        *(
            float(np.max(np.abs(audio.mono().samples)))
            for pair in pairs
            for audio in (pair.reference, pair.synthesis)
        ),
        *(float(np.max(np.abs(clip.audio.mono().samples))) for clip in auditions),
        1.0e-9,
    )
    return float(min(requested, 0.9 / peak))


def _aligned_pair(pair: ReportPair) -> tuple[AudioBuffer, AudioBuffer]:
    reference = resample_audio(pair.reference.mono(), 48000)
    synthesis = resample_audio(pair.synthesis.mono(), 48000)
    length = min(reference.samples.size, synthesis.samples.size)
    first = reference.samples[:length]
    second = synthesis.samples[:length]
    alignment = measure_alignment(first, second, 48000)
    shifted = shift_with_zeros(second, alignment.candidate_lag_samples)
    return AudioBuffer(first, 48000), AudioBuffer(shifted, 48000)


def _onset_pair(
    reference: AudioBuffer, synthesis: AudioBuffer, seconds: float
) -> tuple[AudioBuffer, AudioBuffer]:
    """Extract equal-length pair segments beginning at the reference impact."""
    count = round(seconds * reference.sample_rate)
    onset = detect_impact_onset(reference.samples, reference.sample_rate)
    return (
        AudioBuffer(
            _fixed_segment(reference.samples, onset, count), reference.sample_rate
        ),
        AudioBuffer(
            _fixed_segment(synthesis.samples, onset, count), synthesis.sample_rate
        ),
    )


def _fixed_segment(samples: np.ndarray, start: int, count: int) -> np.ndarray:
    result = np.zeros(count, dtype=np.float64)
    available = max(0, min(count, samples.size - start))
    if available:
        result[:available] = samples[start : start + available]
    return result


def _repeat_with_gap(
    audio: AudioBuffer, gain: float, repetitions: int = 4, gap_seconds: float = 0.150
) -> AudioBuffer:
    gap = round(gap_seconds * audio.sample_rate)
    stride = audio.samples.size + gap
    result = np.zeros(repetitions * stride - gap, dtype=np.float64)
    for index in range(repetitions):
        start = index * stride
        result[start : start + audio.samples.size] = gain * audio.samples
    return AudioBuffer(result, audio.sample_rate)


def _comparison_figure(
    reference: AudioBuffer, synthesis: AudioBuffer, label: str
) -> go.Figure:
    seconds = np.arange(reference.samples.size) / reference.sample_rate
    stride = max(1, reference.samples.size // 12000)
    reference_stft = stft(
        reference.samples, reference.sample_rate, StftConfig(2048, 256, 4096)
    )
    synthesis_stft = stft(
        synthesis.samples, synthesis.sample_rate, StftConfig(2048, 256, 4096)
    )
    reference_db, synthesis_db = _common_magnitude_db(
        reference_stft.magnitude, synthesis_stft.magnitude, -100.0
    )
    attack_reference, attack_synthesis = _onset_pair(reference, synthesis, 0.030)
    attack_config = StftConfig(128, 16, 256)
    reference_attack = stft(
        attack_reference.samples, attack_reference.sample_rate, attack_config
    )
    synthesis_attack = stft(
        attack_synthesis.samples, attack_synthesis.sample_rate, attack_config
    )
    reference_attack_db, synthesis_attack_db = _common_magnitude_db(
        reference_attack.magnitude, synthesis_attack.magnitude, -80.0
    )
    initial_reference, initial_synthesis = _onset_pair(reference, synthesis, 0.100)
    initial_config = StftConfig(256, 32, 512)
    reference_initial = stft(
        initial_reference.samples, initial_reference.sample_rate, initial_config
    )
    synthesis_initial = stft(
        initial_synthesis.samples, initial_synthesis.sample_rate, initial_config
    )
    reference_initial_db, synthesis_initial_db = _common_magnitude_db(
        reference_initial.magnitude, synthesis_initial.magnitude, -80.0
    )
    figure = make_subplots(
        rows=9,
        cols=1,
        shared_xaxes=False,
        subplot_titles=(
            "Aligned waveform",
            "Reference attack, first 30 ms — 128-sample STFT",
            "TriggerFish attack, first 30 ms — 128-sample STFT",
            "Reference initial decay, first 100 ms — 256-sample STFT",
            "TriggerFish initial decay, first 100 ms — 256-sample STFT",
            "Reference spectrogram — shared dB reference",
            "TriggerFish spectrogram — shared dB reference",
            "Signed spectrogram difference — TriggerFish minus reference",
            "Broad-band decay trajectories — solid reference, dashed TriggerFish",
        ),
        vertical_spacing=0.032,
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
    _add_spectrogram(figure, 2, reference_attack, reference_attack_db, -80.0)
    _add_spectrogram(figure, 3, synthesis_attack, synthesis_attack_db, -80.0)
    _add_spectrogram(figure, 4, reference_initial, reference_initial_db, -80.0)
    _add_spectrogram(figure, 5, synthesis_initial, synthesis_initial_db, -80.0)
    _add_spectrogram(figure, 6, reference_stft, reference_db, -100.0)
    _add_spectrogram(figure, 7, synthesis_stft, synthesis_db, -100.0)
    selected = reference_stft.frequencies_hz >= 80.0
    frequencies = reference_stft.frequencies_hz[selected]
    frequency_stride = max(1, int(np.ceil(frequencies.size / 256)))
    time_stride = max(1, int(np.ceil(reference_stft.times_seconds.size / 400)))
    difference = synthesis_db[selected] - reference_db[selected]
    figure.add_trace(
        go.Heatmap(
            x=reference_stft.times_seconds[::time_stride],
            y=frequencies[::frequency_stride],
            z=difference[::frequency_stride, ::time_stride],
            zmin=-30,
            zmax=30,
            zmid=0,
            colorscale="RdBu",
            colorbar={"title": "Δ dB"},
        ),
        row=8,
        col=1,
    )
    figure.update_yaxes(type="log", range=[np.log10(80), np.log10(22000)], row=8, col=1)
    _add_band_decay_traces(figure, reference_stft, synthesis_stft, 9)
    for boundary in (0.015, 0.25, 1.5, 4.0, 8.0):
        figure.add_vline(x=boundary, line_width=1, line_dash="dot", row=9, col=1)
    figure.update_yaxes(range=[-100, 3], title_text="dB", row=9, col=1)
    figure.update_layout(title=label, height=2850, hovermode="x unified")
    figure.update_xaxes(title_text="Seconds", row=1, col=1)
    figure.update_xaxes(title_text="Seconds", row=9, col=1)
    return figure


def _add_spectrogram(
    figure: go.Figure, row: int, result, magnitude_db: np.ndarray, floor_db: float
) -> None:
    selected = result.frequencies_hz >= 80.0
    frequencies = result.frequencies_hz[selected]
    frequency_stride = max(1, int(np.ceil(frequencies.size / 256)))
    time_stride = max(1, int(np.ceil(result.times_seconds.size / 400)))
    figure.add_trace(
        go.Heatmap(
            x=result.times_seconds[::time_stride],
            y=frequencies[::frequency_stride],
            z=magnitude_db[selected][::frequency_stride, ::time_stride],
            zmin=floor_db,
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


def _add_band_decay_traces(figure: go.Figure, reference, synthesis, row: int) -> None:
    bands = (
        ("100–300 Hz", 100.0, 300.0),
        ("300 Hz–1 kHz", 300.0, 1000.0),
        ("1–3 kHz", 1000.0, 3000.0),
        ("3–6 kHz", 3000.0, 6000.0),
        ("6–10 kHz", 6000.0, 10000.0),
        ("10–18 kHz", 10000.0, 18000.0),
    )
    reference_power = _band_power(reference, bands)
    synthesis_power = _band_power(synthesis, bands)
    peak = max(
        float(np.max(reference_power)),
        float(np.max(synthesis_power)),
        np.finfo(float).tiny,
    )
    floor = peak * 1.0e-10
    colors = ("#4c78a8", "#f58518", "#e45756", "#72b7b2", "#54a24b", "#b279a2")
    for band, (name, _, _), color in zip(range(len(bands)), bands, colors):
        for result, power, source, dash in (
            (reference, reference_power, "Reference", "solid"),
            (synthesis, synthesis_power, "TriggerFish", "dash"),
        ):
            values_db = 10.0 * np.log10(np.maximum(power[band], floor) / peak)
            figure.add_trace(
                go.Scattergl(
                    x=result.times_seconds,
                    y=values_db,
                    name=f"{name} — {source}",
                    line={"color": color, "dash": dash},
                ),
                row=row,
                col=1,
            )


def _band_power(result, bands) -> np.ndarray:
    values = []
    for _, low_hz, high_hz in bands:
        selected = (result.frequencies_hz >= low_hz) & (result.frequencies_hz < high_hz)
        values.append(np.sum(result.power[selected], axis=0))
    return np.asarray(values)


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
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:.35rem .5rem;border-bottom:1px solid #444;text-align:right}}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){{text-align:left}}
</style></head><body><h1>{escape(title)}</h1>
<p>Every comparison is real reference versus the current C++ model. There are no old-model renders.</p>
{sections}</body></html>"""
