"""Check the user's actual strong strike, separately from the reference gesture."""

import json
import os
from copy import deepcopy
from pathlib import Path
import numpy as np
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from refine_user_crash_decay import ProtectedDecay


def main():
    root = Path("build/crash-user-decay")
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        saved = verify_candidate(r, root / "coupled")
        original = json.loads(
            (root / "before" / "user-original.fit.json").read_text(encoding="utf8")
        )
        trial = deepcopy(original)
        for node in trial["instrument"]["nodes"]:
            for key in node["parameters"]:
                node["parameters"][key] = saved["parameters"][key]
        audio = [
            r.decode(r.request(command="renderSnapshot", fit=fit, seconds=6)["pcm"])
            for fit in (original, trial)
        ]
        objective = ProtectedDecay(aligned_reference(r, 6), audio[0], r.sample_rate)
        report = dict(
            event=original["controls"]["event"],
            comparison_note="User gesture differs from reference; reference-relative decay is diagnostic only",
            original=objective.diagnostics(audio[0]),
            candidate=objective.diagnostics(audio[1]),
            peak_db=[float(20 * np.log10(max(1e-15, abs(a).max()))) for a in audio],
        )
        write_wav(
            root / "coupled" / "user-gesture.wav", AudioBuffer(audio[1], r.sample_rate)
        )
        (root / "coupled" / "user-gesture.fit.json").write_text(
            json.dumps(trial, indent=2)
        )
        (root / "coupled" / "user-gesture-audit.json").write_text(
            json.dumps(report, indent=2)
        )
        print(json.dumps(report), flush=True)
    finally:
        r.close()


if __name__ == "__main__":
    main()
