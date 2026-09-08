"""Compare one unchanged crash patch against the available edge-velocity grid."""

import argparse
import json
import os
from pathlib import Path

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from audit_metal_strength import measure


def run(args):
    parameters = json.loads(args.source.read_text())["parameters"]
    rows = []
    for velocity in (24, 48, 72, 96, 120):
        renderer = WorkbenchRenderer(
            os.environ["EMSDK_NODE"],
            "crash-standard",
            Path.cwd(),
            cell=dict(articulation="edge", velocity=velocity, repeat=1),
        )
        try:
            directory = args.output / f"v{velocity:03}"
            directory.mkdir(parents=True, exist_ok=True)
            reference = aligned_reference(renderer, 6)
            candidate = renderer.render(parameters, 6)
            for name, samples in (("reference", reference), ("candidate", candidate)):
                write_wav(
                    directory / (name + ".wav"),
                    AudioBuffer(samples, renderer.sample_rate),
                )
            row = dict(
                velocity=velocity,
                event=renderer.metadata["event"],
                reference_info=renderer.metadata["reference"],
                reference=measure(reference, renderer.sample_rate),
                candidate=measure(candidate, renderer.sample_rate),
            )
            rows.append(row)
            print(json.dumps(row), flush=True)
        finally:
            renderer.close()
    (args.output / "velocity-audit.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
