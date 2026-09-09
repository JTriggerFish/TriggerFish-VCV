"""Archive a user fit exactly, plus a checkpoint at the reference's standard strike."""

import argparse
import json
import os
from pathlib import Path

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def run(args):
    saved = json.loads(args.fit.read_text(encoding="utf8"))
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        if saved["reference"]["sha256"] != renderer.metadata["reference"]["sha256"]:
            raise ValueError("Fit reference is not this standard target")
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "user-original.fit.json").write_text(
            json.dumps(saved, indent=2), encoding="utf8"
        )
        rendered = renderer.request(command="renderSnapshot", fit=saved, seconds=6)
        audio = renderer.decode(rendered["pcm"])
        # The bridge performs explicit saved-fit conversion. Do not reintroduce
        # removed controls by flattening the unconverted archived source.
        current = rendered["fit"]
        write_wav(
            args.output / "user-original.wav", AudioBuffer(audio, renderer.sample_rate)
        )
        parameters = dict(renderer.initial)
        parameters.update(
            {
                k: v
                for node in current["instrument"]["nodes"]
                for k, v in node["parameters"].items()
            }
        )
        reference = aligned_reference(renderer, 6)
        checkpoint(
            renderer,
            SpectralBloomLoss(reference, renderer.sample_rate),
            args.output,
            saved["name"] + " — standard strike",
            parameters,
            reference,
            [
                dict(
                    stage="user starting point",
                    source=str(args.fit),
                    user_event=saved["controls"]["event"],
                    standard_event=renderer.metadata["event"],
                )
            ],
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("fit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
