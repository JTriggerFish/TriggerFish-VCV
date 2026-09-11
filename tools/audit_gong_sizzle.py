"""Read-only upper-layer audit: brightness, modulation and loss blind spots.

Writes diagnostic artifacts only, never presets or audition levels. Normalized
envelopes describe texture, not audio gain or overall perceptual quality.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.fft import fft, ifft
from scipy.signal import resample_poly, stft, welch
from scipy.ndimage import gaussian_filter1d

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from refine_gong_layered_bloom import LayeredBloomLoss


def brightness(audio, rate, start=0.5, end=1.5):
    f, power = welch(audio[round(start * rate) : round(end * rate)], rate, nperseg=4096)
    band = (f >= 3000) & (f <= 15000)
    total = max(power[band].sum(), 1e-30)
    return (
        dict(
            centroid_hz=float((f[band] * power[band]).sum() / total),
            above_7k_fraction=float(power[(f >= 7000) & (f <= 15000)].sum() / total),
        ),
        f,
        power,
    )


def blind_spots(rate):
    time = np.arange(6 * rate) / rate
    tone = lambda hz: 0.1 * np.sin(2 * np.pi * hz * time)
    plain = tone(10000)
    modulated = (
        plain * (1 + 0.8 * np.sin(2 * np.pi * 40 * time)) / np.sqrt(1 + 0.8**2 / 2)
    )
    loss = LayeredBloomLoss(plain, plain, rate)

    # Ignore the boundary transient; compare the same broad high-band cells.
    def difference(a, b):
        error = (loss.envelopes(a) - loss.envelopes(b))[4, 2:-1]
        return float(np.sqrt(np.mean(error * error)))

    texture = ModalTextureLoss(plain, rate, centres=[10000])
    return dict(
        within_band_pitch_shift_8_to_13khz_broad_envelope_error_db=difference(
            tone(8000), tone(13000)
        ),
        forty_hz_am_broad_envelope_error_db=difference(plain, modulated),
        forty_hz_am_texture_distance=texture.score(modulated),
        interpretation="Counterexamples to broad-power sufficiency, not a calibrated hearing threshold",
    )


def local_texture(audio, texture):
    """Separate fast 12-kHz fluctuations from its still-mismatched bloom."""
    env = np.abs(ifft(fft(audio, texture.size) * texture.masks[-1])[: len(audio)])
    env = resample_poly(env, 1, texture.step)
    trend = gaussian_filter1d(env, 0.1 * texture.envelope_rate)
    relative = env / np.maximum(trend, 1e-15)
    segment = relative[
        round(0.6 * texture.envelope_rate) : round(1.6 * texture.envelope_rate)
    ]
    f, p = welch(segment, texture.envelope_rate, nperseg=256)
    power = [
        float(p[(f >= lo) & (f < hi)].sum() * (f[1] - f[0]))
        for lo, hi in ((8, 32), (32, 128))
    ]
    return relative, power


def ridge_contrast(audio, rate):
    """Fine spectral contrast after removing broad colour; diagnostic, not loss.

    Average 40 ms of power to reduce FFT speckle, then subtract an 80 Hz-smoothed
    log spectrum. This distinguishes fine ridges from broad EQ without tracking
    their exact frequencies. Fixed analysis windows apply to both signals.
    """
    f, t, z = stft(audio, rate, nperseg=8192, noverlap=8192 - round(0.01 * rate))
    power = gaussian_filter1d(abs(z) ** 2, 2, axis=1)
    db = 10 * np.log10(np.maximum(power, 1e-20))
    detail = db - gaussian_filter1d(db, 80 / (f[1] - f[0]), axis=0)
    region = (t >= 0.6) & (t < 1.6)
    return [
        float(np.sqrt(np.mean(detail[(f >= lo) & (f < hi)][:, region] ** 2)))
        for lo, hi in ((3000, 7000), (7000, 14000))
    ]


def plot(reference, signals, rate, texture, diagnostics):
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            "Upper spectrum, 0.5–1.5 s; absolute levels",
            "Upper spectral centre (3–15 kHz)",
            f"Model − reference modulation power; {len(signals)}-seed mean",
            "12 kHz envelope; slow trend divided out (diagnosis only)",
        ],
    )
    for name, audio, colour in [
        ("Reference", reference, "#eebc59"),
        ("Model", signals[0], "#68b5ed"),
    ]:
        _, f, p = brightness(audio, rate)
        selected = (f >= 3000) & (f <= 15000)
        fig.add_trace(
            go.Scatter(
                x=f[selected],
                y=10 * np.log10(np.maximum(p[selected], 1e-20)),
                name=name,
                line_color=colour,
            ),
            row=1,
            col=1,
        )
        starts = np.arange(0.1, 2.7, 0.1)
        centres = [
            brightness(audio, rate, t, t + 0.3)[0]["centroid_hz"] for t in starts
        ]
        fig.add_trace(
            go.Scatter(
                x=starts + 0.15,
                y=centres,
                name=name,
                showlegend=False,
                line_color=colour,
            ),
            row=1,
            col=2,
        )
        env, _ = local_texture(audio, texture)
        t = np.arange(len(env)) / texture.envelope_rate
        select = (t >= 0.7) & (t <= 1)
        fig.add_trace(
            go.Scatter(
                x=t[select],
                y=env[select],
                name=name,
                showlegend=False,
                line_color=colour,
            ),
            row=2,
            col=2,
        )
    target = texture.target.reshape(-1, 3, 4)
    mean = np.mean(
        [np.array(d["candidate"]).reshape(-1, 3, 4) for d in diagnostics], axis=0
    )
    select = texture.centres >= 3000
    fig.add_trace(
        go.Heatmap(
            x=["2–8 Hz", "8–32 Hz", "32–128 Hz"],
            y=texture.centres[select],
            z=10 * (mean - target)[select, 1, :3],
            zmin=-12,
            zmax=12,
            colorscale="RdBu_r",
            colorbar=dict(title="dB"),
        ),
        row=2,
        col=1,
    )
    fig.update_xaxes(title_text="Hz", row=1, col=1)
    fig.update_yaxes(title_text="dB/Hz", row=1, col=1)
    fig.update_xaxes(title_text="seconds", row=1, col=2)
    fig.update_yaxes(title_text="Hz", row=1, col=2)
    fig.update_yaxes(title_text="carrier band, Hz", row=2, col=1)
    fig.update_xaxes(title_text="seconds", row=2, col=2)
    fig.update_layout(
        template="plotly_dark",
        width=1450,
        height=1050,
        title="Gong sizzle audit — brightness and fine fluctuations, not just bloom loudness",
    )
    return fig


def run(args):
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(r, args.directory)
        ref = read_wav(args.directory / "reference.wav").mono().samples
        signals = [r.render(saved["parameters"], 6, s) for s in args.seeds]
        texture = ModalTextureLoss(ref, r.sample_rate)
        diagnostics = [texture.diagnostics(a) for a in signals]
        report = dict(
            reference_brightness=brightness(ref, r.sample_rate)[0],
            candidate_brightness=[brightness(a, r.sample_rate)[0] for a in signals],
            texture=diagnostics,
            counterexamples=blind_spots(r.sample_rate),
            seeds=args.seeds,
            ridge_contrast_db=dict(
                bands_hz=[[3000, 7000], [7000, 14000]],
                reference=ridge_contrast(ref, r.sample_rate),
                candidate=[ridge_contrast(a, r.sample_rate) for a in signals],
            ),
            preset_modified=False,
            detrended_12k=dict(
                trend_sigma_seconds=0.1,
                region_seconds=[0.6, 1.6],
                modulation_bands_hz=[[8, 32], [32, 128]],
                reference_power=local_texture(ref, texture)[1],
                candidate_power=[local_texture(a, texture)[1] for a in signals],
            ),
        )
        (args.directory / "sizzle-audit.json").write_text(json.dumps(report, indent=2))
        (args.directory / "sizzle.plotly.json").write_text(
            plot(ref, signals, r.sample_rate, texture, diagnostics).to_json()
        )
        print(json.dumps({k: v for k, v in report.items() if k != "texture"}, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--seeds", type=int, nargs="+", default=[1675, 1982, 2586, 3276]
    )
    run(parser.parse_args())
