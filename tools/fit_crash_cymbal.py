"""Fit the production crash engine and emit one reference-vs-model report."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path

from triggerfish_percussion.audio_io import AudioBuffer, read_wav, write_wav
from triggerfish_percussion.crash_fitting import (
    CrashFitCell,
    fit_body,
    fit_level_curve,
)
from triggerfish_percussion.crash_corpus import load_fit_cells
from triggerfish_percussion.crash_model import (
    CrashEvent,
    CrashFit,
    render_crash,
    render_crash_sequence,
)
from triggerfish_percussion.report import (
    AuditionClip,
    ReportPair,
    write_comparison_report,
)


def bitwig_cells(root: Path) -> tuple[CrashFitCell, ...]:
    strengths = (0.25, 0.5, 0.75, 1.0)
    return tuple(
        CrashFitCell(
            f"Bitwig A Custom 18 Medium edge {index:02d}",
            read_wav(root / f"Crash A Custom 18 med {index:02d}.wav", "mean"),
            strength,
            seed=100 + index,
        )
        for index, strength in enumerate(strengths, 1)
    )


def default_bitwig_root() -> Path:
    local = Path(os.environ["LOCALAPPDATA"])
    return local / (
        "Bitwig Studio/installed-packages/5.0/Bitwig/"
        "Acoustic Drums and Percussion/Acoustic Drums/Cymbals/Crashes"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, default=default_bitwig_root())
    parser.add_argument(
        "--cells-manifest",
        type=Path,
        help="isolated SD3 cells.json; omit for the preliminary Bitwig edge fit",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/cymbal-calibration/crash-v1"),
    )
    parser.add_argument("--maximum-evaluations", type=int, default=800)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--load-fit",
        type=Path,
        help="reuse parameters from fit.json and regenerate renders/report",
    )
    arguments = parser.parse_args()
    all_cells = (
        load_fit_cells(arguments.cells_manifest)
        if arguments.cells_manifest
        else bitwig_cells(arguments.reference_root)
    )
    cells = (
        load_fit_cells(arguments.cells_manifest, "fit")
        if arguments.cells_manifest
        else all_cells
    )
    if len(cells) < 2:
        parser.error("the selected corpus needs at least two fitting cells")
    body_cells = _body_anchors(cells)
    if arguments.load_fit:
        loaded = json.loads(arguments.load_fit.read_text(encoding="utf-8"))
        fitted = CrashFit(**loaded["parameters"])
        diagnostics = loaded.get("diagnostics", {}) | {
            "report_regenerated_from": str(arguments.load_fit)
        }
    else:
        initial = fit_level_curve(cells, CrashFit())
        fitted, diagnostics = fit_body(
            body_cells, initial, arguments.maximum_evaluations, arguments.seed
        )
        fitted = fit_level_curve(cells, fitted)
    arguments.output.mkdir(parents=True, exist_ok=True)
    fit_path = arguments.output / "fit.json"
    fit_path.write_text(
        json.dumps(
            {
                "schema": 1,
                "training_role": (
                    "primary SD3 multi-location fit"
                    if arguments.cells_manifest
                    else "preliminary Bitwig edge-velocity fit; SD3 audio pending capture"
                ),
                "cell_count": len(cells),
                "body_anchor_labels": [cell.label for cell in body_cells],
                "parameters": asdict(fitted),
                "diagnostics": diagnostics,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    pairs = []
    report_cells = _report_selection(all_cells)
    for index, cell in enumerate(report_cells, 1):
        samples = render_crash(
            fitted,
            12.0,
            48000,
            cell.strength,
            cell.location,
            cell.hardness,
            cell.seed,
        )
        synthesis = AudioBuffer(samples, 48000)
        write_wav(arguments.output / f"synth-{index:02d}.wav", synthesis)
        pairs.append(ReportPair(cell.label, cell.reference, synthesis))
    report = write_comparison_report(
        tuple(pairs),
        arguments.output / "report.html",
        "Crash cymbal calibration",
        _control_sweeps(fitted),
    )
    print(f"wrote {fit_path}")
    print(f"wrote {report}")
    print(json.dumps(diagnostics, indent=2))


def _body_anchors(cells: tuple[CrashFitCell, ...]) -> tuple[CrashFitCell, ...]:
    """Keep global search tractable while spanning strength and location."""
    groups: dict[tuple[float, float], list[CrashFitCell]] = {}
    for cell in cells:
        groups.setdefault((cell.location, cell.hardness), []).append(cell)
    anchors = []
    for group in groups.values():
        ordered = sorted(group, key=lambda cell: cell.strength)
        anchors.append(ordered[-1])
        if ordered[-1].location == 1.0 and ordered[-1].hardness < 0.8:
            anchors.append(ordered[len(ordered) // 2])
    return tuple(anchors)


def _report_selection(cells: tuple[CrashFitCell, ...]) -> tuple[CrashFitCell, ...]:
    """Show low/middle/high dynamics per control point without an 80-plot page."""
    if len(cells) <= 16:
        return cells
    groups: dict[tuple[float, float], list[CrashFitCell]] = {}
    for cell in cells:
        groups.setdefault((cell.location, cell.hardness), []).append(cell)
    selected = []
    for group in groups.values():
        ordered = sorted(group, key=lambda cell: cell.strength)
        for index in sorted({0, len(ordered) // 2, len(ordered) - 1}):
            selected.append(ordered[index])
    return tuple(selected)


def _control_sweeps(fit: CrashFit) -> tuple[AuditionClip, ...]:
    """Render orthogonal one-second control sweeps through one stateful body."""
    values = (0.15, 0.32, 0.49, 0.66, 0.83, 1.0)
    specifications = (
        ("Velocity, quiet to loud", values, (0.75,) * 6, (0.60,) * 6),
        ("Location, bell to edge", (0.72,) * 6, values, (0.60,) * 6),
        ("Hardness, soft to hard", (0.72,) * 6, (0.75,) * 6, values),
    )
    clips = []
    for clip_index, (label, strengths, locations, hardnesses) in enumerate(
        specifications
    ):
        events = tuple(
            CrashEvent(
                0.25 + event_index,
                strengths[event_index],
                locations[event_index],
                hardnesses[event_index],
                0x53574550 + 16 * clip_index + event_index,
            )
            for event_index in range(6)
        )
        clips.append(
            AuditionClip(
                label,
                AudioBuffer(render_crash_sequence(fit, 10.0, events), 48000),
            )
        )
    return tuple(clips)


if __name__ == "__main__":
    main()
