"""Correct measured ringing through visible controls, preserving output gain."""

import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.ridge_balance_loss import RidgeBalanceLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def main():
    root = Path(__file__).resolve().parents[1]
    direct = os.environ.get("TF_KICK_DIRECT_START") == "1"
    prune = os.environ.get("TF_KICK_PRUNE") == "1"
    noise_only = os.environ.get("TF_KICK_NOISE_OBSERVATION") == "1"
    lowpass = os.environ.get("TF_KICK_LOWPASS") == "1"
    measured = os.environ.get("TF_KICK_MEASURED_BASS") == "1"
    variant = os.environ.get("TF_KICK_VARIANT", "")
    pulse_body = os.environ.get("TF_KICK_PULSE_BODY") == "1"
    if variant not in ("", "short", "rounded"):
        raise ValueError("Unknown controlled kick variation")
    output = root / (
        "build/kick-ridge-pulsebody"
        if pulse_body
        else (
            f"build/kick-ridge-{variant}"
            if variant
            else (
                "build/kick-ridge-measured"
                if measured
                else (
                    "build/kick-ridge-lowpass"
                    if lowpass
                    else (
                        "build/kick-ridge-noise"
                        if noise_only
                        else (
                            "build/kick-ridge-pruned"
                            if prune
                            else (
                                "build/kick-ridge-direct"
                                if direct
                                else "build/kick-ridge-fullband"
                            )
                        )
                    )
                )
            )
        )
    )
    output.mkdir(exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset:]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        baseline = renderer.render(renderer.initial, 1.2)
        loss = RidgeBalanceLoss(reference, baseline, rate)
        search = Search(
            renderer, loss, output, 1.2, "Kick ridge correction", (1449, 1450)
        )
        if prune:
            p = search.parameters
            floor = 10 ** (-72 / 20)
            weights = [
                max(0, 10 ** (p[f"resonance_level_{i}"] / 20) - floor)
                for i in range(16)
            ]
            p["resonance_level"] *= (sum(weights) - weights[4] - weights[5]) / sum(
                weights
            )
            p["resonance_level_4"] = p["resonance_level_5"] = -72
            p.update(
                contact_level=0.1,
                contact_noise_decay_seconds=0.06,
                thump_decay_shape=0.5,
                thump_hold_seconds=0,
            )
        if direct:
            prior = json.loads(
                (root / "build/kick-ridge-fullband/search.json").read_text()
            )
            search.parameters.update(prior["parameters"])
            # Explicit alternate source regime, not removal of contact noise:
            # short broadband direct attack, sustained body supplied by modes.
            search.parameters.update(
                contact_level=0.2,
                contact_noise_level=3,
                contact_noise_decay_seconds=0.035,
                contact_width_seconds=0.003,
            )
        if noise_only:
            prior = json.loads(
                (root / "build/kick-ridge-pruned/search.json").read_text()
            )
            search.parameters.update(prior["parameters"])
            search.parameters.update(
                contact_observation=1,
                contact_level=0.2,
                contact_noise_decay_seconds=0.05,
            )
        if lowpass:
            search.parameters.update(
                equalizer_mode=1, low_cut_hz=5, high_cut_hz=1200, colour_gain_db=0
            )
        if measured:
            prior = json.loads(
                (root / "build/kick-ridge-lowpass/search.json").read_text()
            )
            search.parameters.update(prior["parameters"])
            search.parameters.update(resonance_frequency_7=118, resonance_level_7=-28)
        if variant:
            prior = json.loads(
                (root / "build/kick-ridge-measured/search.json").read_text()
            )
            search.parameters.update(prior["parameters"])
            if variant == "short":
                search.parameters["contact_width_seconds"] = 0.001
            else:
                search.parameters.update(
                    thump_decay_shape=1, thump_hold_seconds=0, thump_decay_seconds=0.24
                )
        if pulse_body:
            prior = json.loads(
                (root / "build/kick-ridge-measured/search.json").read_text()
            )
            search.parameters.update(prior["parameters"])
            search.parameters.update(
                contact_body_drive=1,
                contact_observation=1,
                contact_width_seconds=0.001,
                contact_noise_level=1,
                contact_level=1,
                contact_noise_decay_seconds=0.08,
                resonance_decay_seconds=0.45,
                resonance_decay_tilt=0,
            )
        bounds = dict(
            contact_level=(0, 4),
            contact_noise_level=(0.1, 4),
            contact_noise_decay_seconds=(0.01, 0.25),
            contact_colour=(0, 1),
            thump_level=(1, 4),
            thump_decay_seconds=(0.12, 0.36),
            thump_hold_seconds=(0, 0.03),
            resonance_level=(0.3, 12),
            resonance_decay_seconds=(0.15, 0.8),
            resonance_decay_tilt=(-1, 1),
            thump_pitch_hz=(23, 34),
            thump_pitch_drop_octaves=(0.7, 2.3),
            thump_pitch_fall_seconds=(0.02, 0.09),
            contact_width_seconds=(0.002, 0.025),
        )
        # Modal levels are relative. Freeze one anchor to remove their common
        # scale nullspace; audible bank level remains an explicit fitted control.
        bounds.update({f"resonance_level_{i}": (-60, -4) for i in range(8) if i != 1})
        if prune:
            for i in (4, 5):
                del bounds[f"resonance_level_{i}"]
            bounds["thump_decay_shape"] = (0, 1)
        if lowpass:
            bounds["high_cut_hz"] = (500, 4000)
        if measured:
            bounds.update(
                resonance_frequency_0=(48, 65),
                resonance_frequency_1=(78, 99),
                resonance_frequency_2=(24, 36),
                resonance_frequency_7=(108, 125),
            )
        if variant == "short":
            bounds["contact_width_seconds"] = (0.0002, 0.003)
        if variant == "rounded":
            del bounds["thump_decay_shape"]
            del bounds["thump_hold_seconds"]
        if pulse_body:
            # Noise amount and direct level form a pure product in this routing.
            # Fit its visible observation level; leave source noise amount at 1.
            del bounds["contact_noise_level"]
        search.stage(
            "Resolve spectral excesses and source balance",
            bounds,
            35 if noise_only else 55,
        )
        # Do not force the old frequencies to survive in the fitted answer.
        for i in (3, 4, 5, 6, 7):
            if prune and i in (4, 5):
                continue
            if measured and i == 7:
                continue
            frequency = renderer.initial[f"resonance_frequency_{i}"]
            bounds[f"resonance_frequency_{i}"] = (frequency * 0.85, frequency * 1.15)
        if not variant:
            search.stage(
                "Refine body frequencies with ridge constraints",
                bounds,
                30 if noise_only else 55,
            )
        checks = []
        for seed in (1449, 1450, 1451, 1452):
            checks.append(
                dict(
                    seed=seed,
                    before=loss.diagnostics(
                        renderer.render(renderer.initial, 1.2, seed)
                    ),
                    after=loss.diagnostics(
                        renderer.render(search.parameters, 1.2, seed)
                    ),
                )
            )
        (output / "checks.json").write_text(json.dumps(checks, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
