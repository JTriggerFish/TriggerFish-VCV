"""Local body-envelope experiment; exact Wasm, fixed gain and contact controls.

Candidates are diagnostic artifacts, never automatically published as presets.
The logged finite differences establish which controls influence the objective.
"""

import json
import os
from pathlib import Path

import numpy as np

from triggerfish_percussion.drum_balance_loss import DrumBalanceLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def main():
    root = Path(__file__).resolve().parents[1]
    hold = float(os.environ.get("TF_KICK_HOLD", "0"))
    output = root / f"build/kick-body-hold-{round(hold * 1000)}ms"
    output.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset:]
        reference = np.pad(reference, (0, max(0, round(1.2 * rate) - len(reference))))
        loss = DrumBalanceLoss(reference, rate)
        search = Search(renderer, loss, output, 1.2, "Kick body", (1449, 1450))
        search.parameters["thump_hold_seconds"] = hold
        before = loss.diagnostics(search.audio(search.parameters))
        bounds = {
            "thump_decay_seconds": (0.08, 0.4),
            "thump_level": (1, 4),
            "resonance_decay_seconds": (0.1, 0.7),
            "resonance_decay_tilt": (-0.5, 1),
            "resonance_level": (1, 8),
            **{f"resonance_level_{i}": (-45, 3) for i in range(6)},
        }
        search.stage(
            "Shared damping and body balance; fixed frequencies/contact", bounds, 35
        )
        sparse = dict(search.parameters)
        # The reference retains 90–350 Hz energy after the two bass sources fade.
        # Test two existing spare handles instead of stretching every mode's T60.
        for i in (6, 7):
            search.parameters[f"resonance_level_{i}"] = -35
            bounds[f"resonance_level_{i}"] = (-60, -8)
        search.stage("Two upper-bass handles; shared damping only", bounds, 35)
        search.screen_candidates(
            "Keep original sparse alternative eligible", [("six modes", sparse)]
        )
        report = dict(
            before=before,
            after=loss.diagnostics(search.audio(search.parameters)),
            parameters=search.parameters,
            metadata=renderer.metadata,
        )
        (output / "comparison.json").write_text(json.dumps(report, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
