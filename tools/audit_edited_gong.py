"""Holdout gestures and seeds for a saved-series refinement, without level matching."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from refine_edited_gong import Comparison
from study_saved_gong_texture import texture_features


def run(args):
    source = args.directory / "original.fit.json"
    if not source.exists():
        source = args.directory / "source.fit.json"
    original = json.loads(source.read_text())
    candidate = json.loads((args.directory / "candidate.fit.json").read_text())
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = SavedFitRenderer(r, original)
        parameters = saved.parameters(candidate)
        changes = {k for k, v in parameters.items() if v != saved.initial[k]}
        if not changes <= {
            "bloom_rate",
            "body_decay_seconds_7",
            "field_motion_depth",
            "field_packet_spread",
        }:
            raise ValueError(
                "The selected refinement changed more than shared dynamics"
            )
        if candidate["controls"] != original["controls"]:
            raise ValueError(
                "The selected refinement changed the saved controls/gesture"
            )
        reference = aligned_reference(r, 6)
        texture = ModalTextureLoss(
            reference, r.sample_rate, centres=[2500, 3500, 4500, 6500, 9500, 12500]
        )
        reference_texture = texture_features(texture, reference)
        rows = []
        for name, event in (
            ("saved strike", original["controls"]["event"]),
            ("reference strike", r.metadata["event"]),
        ):
            for seed in (2673, 3911):
                before = saved.render(saved.initial, 6, seed, event)
                after = saved.render(parameters, 6, seed, event)
                loss = Comparison(reference, before, r.sample_rate)
                row = dict(
                    event=name,
                    seed=seed,
                    before=loss.metrics(before),
                    after=loss.metrics(after),
                )
                for key, audio in (("before", before), ("after", after)):
                    row[key]["texture_db"] = float(
                        np.sqrt(
                            np.mean(
                                (texture_features(texture, audio) - reference_texture)
                                ** 2
                            )
                        )
                    )
                rows.append(row)
                print(json.dumps(row), flush=True)
        energies = []
        for strength in (0.3, 0.5, 0.76, 1):
            audio = saved.render(parameters, 6, event=dict(strength=strength))
            if not np.all(np.isfinite(audio)):
                raise ValueError("Nonfinite rendered output")
            energies.append(
                dict(
                    strength=strength,
                    energy=float(np.sum(audio**2) / r.sample_rate),
                    peak_db=float(20 * np.log10(max(abs(audio).max(), 1e-20))),
                )
            )
        if not all(b["energy"] > a["energy"] for a, b in zip(energies, energies[1:])):
            raise ValueError("Increasing strength did not increase output energy")
        report = dict(holdouts=rows, velocities=energies, changed=sorted(changes))
        (args.directory / "holdout.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(dict(velocities=energies)), flush=True)
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    run(parser.parse_args())
