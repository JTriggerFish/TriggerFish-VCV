"""Split one rendered cymbal sweep into deterministic fitting cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.crash_corpus import isolate_capture, write_cells_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audio", type=Path, required=True)
    parser.add_argument("--sweep-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration-seconds", type=float, default=10.0)
    parser.add_argument("--pre-roll-seconds", type=float, default=0.05)
    arguments = parser.parse_args()
    if arguments.duration_seconds <= 0 or arguments.pre_roll_seconds < 0:
        parser.error("duration must be positive and pre-roll non-negative")
    sweep = json.loads(arguments.sweep_manifest.read_text(encoding="utf-8"))
    cells = isolate_capture(
        read_wav(arguments.audio, "mean"),
        sweep,
        arguments.output,
        arguments.duration_seconds,
        pre_roll_seconds=arguments.pre_roll_seconds,
    )
    manifest = arguments.output / "cells.json"
    write_cells_manifest(manifest, cells, arguments.audio)
    print(f"wrote {len(cells)} isolated cells and {manifest}")


if __name__ == "__main__":
    main()
