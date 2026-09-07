"""Fit low-body envelope without redistributing other modal observations."""

import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from triggerfish_percussion.drum_balance_loss import DrumBalanceLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.ridge_balance_loss import RidgeBalanceLoss
from triggerfish_percussion.workbench_search import Search


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "build/kick-ridge-lowbody"
    output.mkdir(exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        if os.environ.get("TF_KICK_COLOUR_ONLY") == "1":
            finish_colour(renderer, root)
            return
        prior = json.loads(
            (root / "build/kick-ridge-pulsebody/search.json").read_text()
        )
        initial = prior["parameters"]
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset:]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        loss = DrumBalanceLoss(reference, rate)
        keys = [
            "thump_level",
            "thump_decay_seconds",
            "thump_hold_seconds",
            "thump_decay_shape",
            "thump_pitch_hz",
            "resonance_frequency_2",
            "resonance_level_2",
        ]
        low = np.array([1, 0.1, 0, 0, 24, 25, -50.0])
        high = np.array([4, 0.5, 0.04, 1, 29, 29.5, -8.0])
        floor = 10 ** (-72 / 20)
        weight = lambda p: sum(
            max(0, 10 ** (p[f"resonance_level_{i}"] / 20) - floor) for i in range(16)
        )
        original_weight = weight(initial)
        cache = {}

        def unpack(x):
            p = dict(initial, **dict(zip(keys, (low + x * (high - low)).tolist())))
            # Static authoring conversion, exported in the actual patch. Every
            # other mode retains its absolute coefficient despite normalization.
            p["resonance_level"] = (
                initial["resonance_level"] * weight(p) / original_weight
            )
            return p

        def residual_parameters(p):
            cache_key = tuple(p[k] for k in keys)
            if cache_key not in cache:
                audio = renderer.render(p, 1.2, 1449)
                power = loss.envelopes(audio)
                delta = (
                    10 * np.log10(np.maximum(power[:2], loss.floor)) - loss.target[:2]
                )
                selected = loss.times < 0.45
                w = loss.weight[:2, selected]
                bass = (delta[:, selected] * np.sqrt(w / w.sum())).ravel()
                cache[cache_key] = np.concatenate((bass, 0.35 * loss.residual(audio)))
            return cache[cache_key]

        def residual(x):
            return residual_parameters(unpack(x))

        def jacobian(x):
            columns = []
            for i in range(len(x)):
                a, b = x.copy(), x.copy()
                a[i], b[i] = max(0, x[i] - 0.005), min(1, x[i] + 0.005)
                columns.append((residual(b) - residual(a)) / (b[i] - a[i]))
            return np.array(columns).T

        before = float(np.linalg.norm(residual_parameters(initial)))
        start = dict(
            initial,
            thump_level=3,
            thump_decay_seconds=0.26,
            thump_hold_seconds=0,
            resonance_frequency_2=27,
            resonance_level_2=-18,
        )
        x = np.clip((np.array([start[k] for k in keys]) - low) / (high - low), 0, 1)
        result = least_squares(
            residual,
            x,
            bounds=(0, 1),
            jac=jacobian,
            max_nfev=65,
            ftol=0.001,
            xtol=0.001,
            gtol=0.001,
        )
        parameters = unpack(result.x)
        after = float(np.linalg.norm(residual_parameters(parameters)))
        fit = renderer.request(
            command="snapshot", parameters=parameters, name="Kick — low body trial"
        )["fit"]
        (output / "candidate.fit.json").write_text(json.dumps(fit, indent=2))
        (output / "search.json").write_text(
            json.dumps(
                dict(
                    parameters=parameters,
                    before=before,
                    after=after,
                    evaluations=len(cache),
                    keys=keys,
                    bounds=[low.tolist(), high.tolist()],
                    metadata=renderer.metadata,
                ),
                indent=2,
            )
        )
        print(
            json.dumps(dict(before=before, after=after, evaluations=len(cache))),
            flush=True,
        )
    finally:
        renderer.close()


def finish_colour(renderer, root):
    """One broad observation peak and low-pass; no pole or source-envelope fit."""
    output = root / "build/kick-ridge-colour"
    output.mkdir(exist_ok=True)
    parameters = json.loads(
        (root / "build/kick-ridge-lowbody/search.json").read_text()
    )["parameters"]
    rate = renderer.sample_rate
    onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
    reference = renderer.reference[onset:]
    reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
    loss = RidgeBalanceLoss(reference, renderer.render(renderer.initial, 1.2), rate)
    search = Search(renderer, loss, output, 1.2, "Kick — pulse body", (1449, 1450))
    search.parameters = parameters
    search.stage(
        "Broad radiation colour; unchanged modal poles and source envelopes",
        dict(colour_gain_db=(0, 10), high_cut_hz=(500, 2000), contact_level=(0.3, 4)),
        35,
    )
    # The centre has zero influence while gain is zero. Recheck it after the
    # first stage activates colour instead of freezing that direction forever.
    search.stage(
        "Refine the now-audible colour centre",
        dict(
            colour_frequency_hz=(300, 1800),
            colour_gain_db=(0, 10),
            high_cut_hz=(500, 2000),
            contact_level=(0.3, 4),
        ),
        30,
    )


if __name__ == "__main__":
    main()
