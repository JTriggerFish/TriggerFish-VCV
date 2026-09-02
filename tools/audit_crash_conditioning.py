"""Audit 0--100 ms crash parameter leverage before fitting."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

from fit_crash_cymbal import bitwig_cells, default_bitwig_root, _select_fit_cell
from triggerfish_percussion.crash_fit_conditioning import (
    morris_screening_diagnostics,
    parameter_conditioning_diagnostics,
)
from triggerfish_percussion.crash_model import CrashFit
from triggerfish_percussion.crash_fit_parameters import SINGLE_HIT_ATTACK_PARAMETERS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument(
        "--fit-cell",
        required=True,
        help="exact Bitwig cell label or one-based index",
    )
    parser.add_argument("--candidate-stage")
    parser.add_argument("--reference-root", type=Path, default=default_bitwig_root())
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--method", choices=("local", "morris", "both"), default="both")
    parser.add_argument("--trajectories", type=int, default=8)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    payload = json.loads(arguments.fit.read_text(encoding="utf-8"))
    fit = CrashFit(**payload["parameters"])
    if arguments.candidate_stage:
        stages = [
            stage
            for stage in payload.get("diagnostics", {})
            .get("causal_fit", {})
            .get("stages", ())
            if stage.get("stage") == arguments.candidate_stage
            and stage.get("candidate_parameters")
        ]
        if not stages:
            parser.error("requested candidate stage is absent from fit diagnostics")
        fit = replace(fit, **stages[-1]["candidate_parameters"])
    source_cells = bitwig_cells(arguments.reference_root)
    cells = (_select_fit_cell(source_cells, arguments.fit_cell, parser),)
    diagnostics = {}
    if arguments.method in ("local", "both"):
        diagnostics["local"] = parameter_conditioning_diagnostics(
            cells, fit, parameters=SINGLE_HIT_ATTACK_PARAMETERS
        )
    if arguments.method in ("morris", "both"):
        diagnostics["morris"] = morris_screening_diagnostics(
            cells,
            fit,
            parameters=SINGLE_HIT_ATTACK_PARAMETERS,
            trajectories=arguments.trajectories,
            workers=arguments.workers,
        )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
