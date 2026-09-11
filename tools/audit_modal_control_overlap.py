"""Read-only control-reduction experiments using the actual workbench DSP.

No preset writes, gain matching or automatic perceptual acceptance. Compare
ablations with seed variation, and isolated rings with all transfer disabled.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import hilbert, welch

from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def colour(audio, rate):
    """Fixed broad spectral bands, intentionally not a complete texture metric."""
    f, p = welch(audio[round(0.1 * rate) :], rate, nperseg=8192)
    edges = [40, 250, 800, 2500, 5000, 9000, 15000]
    return np.array(
        [
            10 * np.log10(max(1e-30, p[(f >= a) & (f < b)].sum()))
            for a, b in zip(edges[:-1], edges[1:])
        ]
    )


def variants(p):
    yield "baseline", dict(p)
    for name, key in [
        ("no beating", "field_beat_depth"),
        ("no wander", "field_wander_hz"),
        ("no blur", "field_phase_bandwidth"),
        ("no ridge movement", "field_motion_depth"),
    ]:
        if p[key] > 0:
            yield name, dict(p, **{key: 0})
    sharing = p["field_motion_sharing"]
    if p["field_motion_depth"] > 0:
        yield "independent motion, variance matched", dict(
            p,
            field_motion_sharing=0,
            field_motion_depth=p["field_motion_depth"] * np.hypot(1 - sharing, sharing),
        )
    yield "fixed beat slope 0.25", dict(p, field_beat_rate_tilt=0.25)
    if p["field_wander_hz"] > 0 and p["field_motion_depth"] == 0:
        for depth in (0.6, 1.2, 2.4):
            yield f"wander replaced by slow shimmer {depth}", dict(
                p,
                field_wander_hz=0,
                field_motion_depth=depth,
                field_motion_rate=p["field_wander_rate"],
                field_motion_sharing=0,
            )


def preset_study(name, output):
    fit = json.loads(
        Path(f"workbench/web/{name}_calibration.fit.json").read_text(encoding="utf-8")
    )
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], name + "-standard", Path.cwd()
    )
    try:
        saved = SavedFitRenderer(renderer, fit)
        seeds = [fit["controls"]["event"]["seed"], 1982]
        bases = [saved.render(saved.initial, 4, seed) for seed in seeds]
        losses = [
            ModalTextureLoss(
                audio, renderer.sample_rate, centres=[125, 500, 1500, 4000, 8000, 12000]
            )
            for audio in bases
        ]
        base_colours = [colour(audio, renderer.sample_rate) for audio in bases]
        rows = []
        for title, parameters in variants(saved.initial):
            metrics = []
            for i, seed in enumerate(seeds):
                audio = (
                    bases[i]
                    if title == "baseline"
                    else saved.render(parameters, 4, seed)
                )
                metrics.append(
                    dict(
                        band_change_db=(
                            colour(audio, renderer.sample_rate) - base_colours[i]
                        ).tolist(),
                        texture_distance=losses[i].score(audio),
                        pcm_error_rms=float(np.sqrt(np.mean((audio - bases[i]) ** 2))),
                    )
                )
            rows.append(dict(name=title, metrics=metrics))
            print(name, title, json.dumps(metrics), flush=True)
        result = dict(
            fit_id=fit["id"],
            seeds=seeds,
            rows=rows,
            parameters=saved.initial,
            renderer_sha256=renderer.metadata["rendererSha256"],
            seed_variation=dict(
                texture_distance=losses[0].score(bases[1]),
                band_change_db=(base_colours[1] - base_colours[0]).tolist(),
            ),
        )
        (output / f"{name}.json").write_text(json.dumps(result, indent=2))
    finally:
        renderer.close()


def isolated_study(output):
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        fit = json.loads(
            Path("workbench/web/gong_calibration.fit.json").read_text(encoding="utf-8")
        )
        saved = SavedFitRenderer(renderer, fit)
        p = dict(
            saved.initial,
            body_tune=1,
            bloom_rate=0,
            output_eq_enabled=0,
            field_distribution=3,
            field_turbulence=1,
            field_turbulence_slope=0,
            field_satellite_density=0,
            field_beat_depth=0,
            field_doublet_split=2,
            field_beat_rate_tilt=0,
            field_motion_depth=0,
            field_motion_rate=40,
            field_motion_sharing=0,
            field_phase_bandwidth=0,
            field_phase_tilt=0,
            field_wander_hz=0,
            field_wander_rate=1,
            body_decay_seconds_0=30,
            body_decay_seconds_7=30,
            direct_gain=0,
        )
        p.update({f"resolved_level_{i}": -72 for i in range(32)})
        p.update({key: 0 for key in p if key.startswith("body_decay_active_")})
        p.update(resolved_level_0=-6, resolved_frequency_0=1000)
        cases = {
            "plain": {},
            "beat 2 Hz": {"field_beat_depth": 0.3},
            "blur narrow": {"field_phase_bandwidth": 0.02},
            "blur broad": {"field_phase_bandwidth": 0.2},
            "wander slow": {"field_wander_hz": 0.3},
            "wander large": {"field_wander_hz": 3},
            "shimmer gentle": {"field_motion_depth": 0.6},
            "shimmer strong": {"field_motion_depth": 1.5},
            "shimmer maximum": {"field_motion_depth": 3},
            "shimmer slow": {"field_motion_depth": 3, "field_motion_rate": 1},
            "paired, independent": {"field_beat_depth": 0.3, "field_motion_depth": 1.5},
            "paired, together": {
                "field_beat_depth": 0.3,
                "field_motion_depth": 1.5,
                "field_motion_sharing": 1,
            },
            "paired, slow wander": {"field_beat_depth": 0.3, "field_wander_hz": 0.3},
            "paired, slow shimmer": {
                "field_beat_depth": 0.3,
                "field_motion_depth": 1.2,
                "field_motion_rate": 1,
            },
        }
        rows = []
        for name, edits in cases.items():
            audio = saved.render(dict(p, **edits), 8)
            rate = renderer.sample_rate
            t = np.arange(len(audio)) / rate
            # Remove the known decay for this diagnostic, not audition gain.
            segment = (audio * 10 ** (3 * t / 30))[round(rate) : round(7 * rate)]
            f, power = welch(segment, rate, nperseg=round(2 * rate))
            analytic = hilbert(segment)
            core = float(power[abs(f - 1000) <= 1].sum() / power.sum())
            env = abs(analytic)[round(0.1 * rate) : -round(0.1 * rate)]
            ef, ep = welch(env, rate, nperseg=round(4 * rate))
            rows.append(
                dict(
                    name=name,
                    carrier_band_fraction=core,
                    envelope_cv=float(np.std(env) / np.mean(env)),
                    periodic_beat_fraction=float(
                        ep[abs(ef - 2) <= 0.25].sum()
                        / max(1e-30, ep[(ef >= 0.5) & (ef <= 10)].sum())
                    ),
                    wide_wings_fraction=float(
                        power[abs(f - 1000) > 10].sum() / power.sum()
                    ),
                )
            )
        (output / "isolated.json").write_text(json.dumps(rows, indent=2))
        print(json.dumps(rows, indent=2), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=["gong", "crash", "ride", "isolated"])
    args = parser.parse_args()
    output = Path("build/modal-control-overlap")
    output.mkdir(parents=True, exist_ok=True)
    if args.only in (None, "isolated"):
        isolated_study(output)
    for instrument in ("gong", "crash", "ride"):
        if args.only in (None, instrument):
            preset_study(instrument, output)
