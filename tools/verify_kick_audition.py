"""Verify exact reload/retriggers and retain the full-match gate's result."""

import json
import os
from pathlib import Path

import numpy as np
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_rerender import prepare_current_render
from triggerfish_percussion.kick_quality_checks import check_kick_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from kick_playability import check_playability, check_resonance_mix


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "build/kick-final-audition"
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        prepare_current_render(renderer, root / "build/kick-ridge-colour", output)
        saved = json.loads((output / "search.json").read_text())
        count = round(saved["duration_seconds"] * renderer.sample_rate)
        onset = round(
            renderer.metadata["reference"]["cell"]["onset_seconds"]
            * renderer.sample_rate
        )
        reference = renderer.reference[onset : onset + count]
        reference = np.pad(reference, (0, count - len(reference)))
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, output)
        quality = check_kick_candidate(output, renderer=renderer)
        try:
            playability = check_playability(renderer, saved["parameters"], output)
        except RuntimeError:
            for label, overrides in (
                ("bypass", dict(equalizer_mode=0)),
                ("highpass20", dict(low_cut_hz=20)),
                ("highpass5", dict(low_cut_hz=5)),
            ):
                path = output / label
                path.mkdir(exist_ok=True)
                try:
                    check_resonance_mix(
                        renderer, dict(saved["parameters"], **overrides), path
                    )
                    print(label, "affine check passed", flush=True)
                except RuntimeError as error:
                    print(label, str(error), flush=True)
            raise
        result = dict(
            exact_reload=True,
            reference_identity_verified=True,
            full_match_eligible=quality["eligible"],
            listening_approved=False,
            playability=playability,
        )
        (output / "audition-verification.json").write_text(json.dumps(result, indent=2))
        print(json.dumps(result, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
