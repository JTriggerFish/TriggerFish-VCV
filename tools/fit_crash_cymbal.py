"""Fit the production crash engine and emit one reference-vs-model report."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
from time import perf_counter

from triggerfish_percussion.audio_io import AudioBuffer, read_wav, write_wav
from triggerfish_percussion.crash_fitting import (
    CAUSAL_FIT_POLICIES,
    CrashFitCell,
    acceptance_diagnostics,
    fit_causal_model,
    fit_sparse_modes,
    fit_sparse_projection,
    parameter_influence_diagnostics,
)
from triggerfish_percussion.crash_corpus import load_fit_cells
from triggerfish_percussion.crash_fit_parameters import (
    CAUSAL_STAGES,
    SCREENED_ATTACK_STAGES,
    SCREENED_INITIAL_DECAY_STAGES,
    fit_parameter_value,
    replace_fit_parameters,
    single_hit_stages,
)
from triggerfish_percussion.crash_fit_spectral_profile import (
    refine_initial_spectral_profile,
)
from triggerfish_percussion.crash_model import (
    CrashEvent,
    CrashFit,
    render_crash,
    render_crash_sequence,
)
from triggerfish_percussion.report import (
    AuditionClip,
    ReportPair,
    velocity_grid_audition_gain,
    write_comparison_report,
)

FIT_SCHEMA = 5


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
        help="isolated private-corpus cells.json; omit for the preliminary licensed edge fit",
    )
    parser.add_argument(
        "--fit-cell",
        required=True,
        help="exact cell label or one-based cell index; every fit targets one recording",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/cymbal-calibration/crash-v1"),
    )
    parser.add_argument("--maximum-evaluations", type=int, default=800)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument(
        "--fit-policy",
        choices=tuple(CAUSAL_FIT_POLICIES),
        default="strict",
        help="strict causal gates or an explicit whole-100ms attack trade-off",
    )
    parser.add_argument(
        "--maximum-prefix-seconds",
        type=float,
        help="stop the causal bootstrap after this prefix (for focused iteration)",
    )
    parser.add_argument(
        "--start-stage",
        choices=tuple(
            dict.fromkeys(
                stage.name
                for stages in (CAUSAL_STAGES, SCREENED_INITIAL_DECAY_STAGES)
                for stage in stages
            )
        ),
        help="resume at this stage after verifying every earlier absolute gate",
    )
    parser.add_argument(
        "--end-stage",
        choices=tuple(
            dict.fromkeys(
                stage.name
                for stages in (CAUSAL_STAGES, SCREENED_INITIAL_DECAY_STAGES)
                for stage in stages
            )
        ),
        help="stop after this stage for a controlled bootstrap diagnostic",
    )
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--load-fit",
        type=Path,
        help="reuse parameters from fit.json and regenerate renders/report",
    )
    parser.add_argument(
        "--initial-fit",
        type=Path,
        help="seed a new optimization from an earlier fit.json",
    )
    parser.add_argument(
        "--initial-candidate-stage",
        help="overlay the last candidate from this diagnostic stage onto --initial-fit",
    )
    parser.add_argument(
        "--initial-candidate-fit",
        type=Path,
        help="read the diagnostic candidate from a different fit.json",
    )
    parser.add_argument(
        "--initial-candidate-blend",
        type=float,
        default=1.0,
        help="interpolate from --initial-fit to the selected candidate (0..1)",
    )
    parser.add_argument(
        "--screened-attack",
        action="store_true",
        help="fit Morris-screened 100 ms blocks followed by a hard-gated polish",
    )
    parser.add_argument(
        "--screened-initial-decay",
        action="store_true",
        help="fit 100 ms balance, colour, dispersion, and loss blocks before joint polish",
    )
    parser.add_argument(
        "--refine-initial-spectrum",
        action="store_true",
        help="regularize one 100 ms object profile while freezing source balance",
    )
    parser.add_argument(
        "--skip-influence",
        action="store_true",
        help="omit the finite-difference causal influence audit",
    )
    parser.add_argument(
        "--skip-pole-initialization",
        action="store_true",
        help="reuse the seed/default sparse poles for a fast causal-fit iteration",
    )
    parser.add_argument(
        "--recompute-diagnostics",
        action="store_true",
        help="rerun acceptance and influence audits when loading an existing fit",
    )
    arguments = parser.parse_args()
    fit_sources = (
        arguments.load_fit,
        arguments.initial_fit,
    )
    if sum(source is not None for source in fit_sources) > 1:
        parser.error("select only one fit source or refinement mode")
    if arguments.initial_candidate_stage and not arguments.initial_fit:
        parser.error("--initial-candidate-stage requires --initial-fit")
    if arguments.initial_candidate_fit and not arguments.initial_candidate_stage:
        parser.error("--initial-candidate-fit requires --initial-candidate-stage")
    if not 0.0 <= arguments.initial_candidate_blend <= 1.0:
        parser.error("--initial-candidate-blend must lie between zero and one")
    if arguments.screened_attack and arguments.start_stage:
        parser.error("--screened-attack cannot be combined with --start-stage")
    if arguments.screened_attack and arguments.screened_initial_decay:
        parser.error("select only one screened fitting workflow")
    if arguments.refine_initial_spectrum and not arguments.initial_fit:
        parser.error("--refine-initial-spectrum requires --initial-fit")
    if arguments.refine_initial_spectrum and any(
        (
            arguments.screened_attack,
            arguments.screened_initial_decay,
            arguments.start_stage,
            arguments.end_stage,
            arguments.maximum_prefix_seconds,
        )
    ):
        parser.error("spectral-profile refinement is a complete focused workflow")
    all_cells = (
        load_fit_cells(arguments.cells_manifest)
        if arguments.cells_manifest
        else bitwig_cells(arguments.reference_root)
    )
    available_cells = (
        load_fit_cells(arguments.cells_manifest, "fit")
        if arguments.cells_manifest
        else all_cells
    )
    source_cell = _select_fit_cell(available_cells, arguments.fit_cell, parser)
    cells = (source_cell,)
    if arguments.load_fit:
        loaded = json.loads(arguments.load_fit.read_text(encoding="utf-8"))
        _require_loadable_schema(
            loaded, arguments.load_fit, parser, arguments.recompute_diagnostics
        )
        fitted = CrashFit(**loaded["parameters"])
        legacy_schema = loaded.get("schema") != FIT_SCHEMA
        diagnostics = (
            {}
            if legacy_schema and arguments.recompute_diagnostics
            else dict(loaded.get("diagnostics", {}))
        )
        diagnostics["report_regenerated_from"] = str(arguments.load_fit)
        if legacy_schema:
            diagnostics["legacy_loss_audit"] = {
                "source_schema": loaded.get("schema"),
                "status": "rejected; causal diagnostics used obsolete losses",
            }
    else:
        initial_parameters = CrashFit()
        if arguments.initial_fit:
            seed_fit = json.loads(arguments.initial_fit.read_text(encoding="utf-8"))
            _require_current_schema(seed_fit, arguments.initial_fit, parser)
            initial_parameters = CrashFit(**seed_fit["parameters"])
            if arguments.initial_candidate_stage:
                candidate_fit = seed_fit
                if arguments.initial_candidate_fit:
                    candidate_fit = json.loads(
                        arguments.initial_candidate_fit.read_text(encoding="utf-8")
                    )
                    _require_current_schema(
                        candidate_fit, arguments.initial_candidate_fit, parser
                    )
                candidates = [
                    stage
                    for stage in candidate_fit.get("diagnostics", {})
                    .get("causal_fit", {})
                    .get("stages", ())
                    if stage.get("stage") == arguments.initial_candidate_stage
                    and stage.get("candidate_parameters")
                ]
                if not candidates:
                    parser.error(
                        "initial fit has no recorded candidate for "
                        f"{arguments.initial_candidate_stage!r}"
                    )
                selected = min(
                    candidates,
                    key=lambda candidate: float(
                        candidate.get("candidate_composite_objective", float("inf"))
                    ),
                )
                blend = arguments.initial_candidate_blend
                values = {
                    name: fit_parameter_value(initial_parameters, name)
                    + blend * (value - fit_parameter_value(initial_parameters, name))
                    for name, value in selected["candidate_parameters"].items()
                }
                initial_parameters = replace_fit_parameters(initial_parameters, values)
        if arguments.refine_initial_spectrum or arguments.skip_pole_initialization:
            modal_fit = initial_parameters
            modal_diagnostics = {
                "skipped": True,
                "source": "schema-5 seed" if arguments.initial_fit else "defaults",
            }
        else:
            print("extracting persistent pole initialization", flush=True)
            modal_start = perf_counter()
            modal_fit, modal_diagnostics = fit_sparse_modes(cells, initial_parameters)
            modal_fit = fit_sparse_projection(cells, modal_fit)
            modal_diagnostics["elapsed_seconds"] = perf_counter() - modal_start
            print(
                "persistent pole initialization completed in "
                f"{modal_diagnostics['elapsed_seconds']:.1f} s",
                flush=True,
            )
        if arguments.refine_initial_spectrum:
            fitted, causal_diagnostics = refine_initial_spectral_profile(
                source_cell,
                modal_fit,
                arguments.maximum_evaluations,
                progress=lambda message: print(message, flush=True),
            )
        else:
            stages = _selected_stages(arguments, parser)
            fitted, causal_diagnostics = fit_causal_model(
                cells,
                modal_fit,
                arguments.maximum_evaluations,
                arguments.seed,
                stages,
                arguments.start_stage,
                CAUSAL_FIT_POLICIES[arguments.fit_policy],
                progress=lambda message: print(message, flush=True),
                workers=arguments.workers,
            )
        diagnostics = {
            "causal_fit": causal_diagnostics,
            "sparse_pole_initialization": modal_diagnostics,
        }
    reuse_diagnostics = arguments.load_fit and not arguments.recompute_diagnostics
    if not arguments.skip_influence and not reuse_diagnostics:
        diagnostics["parameter_influence"] = parameter_influence_diagnostics(
            max(cells, key=lambda cell: cell.strength), fitted
        )
    if not reuse_diagnostics or "acceptance" not in diagnostics:
        diagnostics["acceptance"] = acceptance_diagnostics(cells, fitted)
    arguments.output.mkdir(parents=True, exist_ok=True)
    fit_path = arguments.output / "fit.json"
    fit_path.write_text(
        json.dumps(
            {
                "schema": FIT_SCHEMA,
                "fit_method": diagnostics.get("causal_fit", {}).get(
                    "method", "cumulative-causal-prefix-v3"
                ),
                "training_role": (
                    "single-hit private-corpus parameter-grid cell"
                    if arguments.cells_manifest
                    else "single-hit preliminary Bitwig parameter-grid cell"
                ),
                "cell_count": len(cells),
                "fit_cell": {
                    "label": source_cell.label,
                    "source_strength_coordinate": source_cell.strength,
                    "location_coordinate": source_cell.location,
                    "hardness_coordinate": source_cell.hardness,
                    "seed": source_cell.seed,
                    "model_strength_during_fit": source_cell.strength,
                },
                "parameters": asdict(fitted),
                "diagnostics": diagnostics,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    pairs = []
    report_cells = tuple(
        cell
        for cell in all_cells
        if cell.strength == source_cell.strength
        and cell.location == source_cell.location
        and cell.hardness == source_cell.hardness
    )
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
        f"Crash cymbal calibration — {source_cell.label}",
        _control_sweeps(fitted),
        diagnostics["acceptance"],
        diagnostics.get("causal_fit"),
        audition_gain=velocity_grid_audition_gain(
            tuple(
                (cell.strength, cell.reference)
                for cell in all_cells
                if cell.location == source_cell.location
                and cell.hardness == source_cell.hardness
            )
        ),
    )
    print(f"wrote {fit_path}")
    print(f"wrote {report}")
    print(json.dumps(diagnostics, indent=2))


def _selected_stages(arguments, parser):
    if arguments.screened_attack:
        stages = SCREENED_ATTACK_STAGES
    elif arguments.screened_initial_decay:
        stages = SCREENED_INITIAL_DECAY_STAGES
    else:
        stages = CAUSAL_STAGES
    stages = single_hit_stages(stages)
    if arguments.maximum_prefix_seconds is not None:
        stages = tuple(
            stage
            for stage in stages
            if stage.end_seconds <= arguments.maximum_prefix_seconds
        )
        if not stages:
            parser.error("maximum prefix is shorter than the first causal stage")
    if arguments.end_stage:
        names = tuple(stage.name for stage in stages)
        if arguments.end_stage not in names:
            parser.error("end stage is outside the selected fitting workflow")
        stages = stages[: names.index(arguments.end_stage) + 1]
    if arguments.start_stage and not any(
        stage.name == arguments.start_stage for stage in stages
    ):
        parser.error("start stage lies beyond the selected maximum prefix")
    return stages


def _require_loadable_schema(
    payload: dict[str, object],
    path: Path,
    parser: argparse.ArgumentParser,
    recompute_diagnostics: bool,
) -> None:
    schema = payload.get("schema")
    if schema == FIT_SCHEMA:
        return
    if schema == 4 and recompute_diagnostics:
        return
    parser.error(
        f"{path} uses crash-fit schema {schema!r}; schema 4 may only be loaded "
        "with --recompute-diagnostics for a read-only audit, and older fits must "
        "be rerun"
    )


def _require_current_schema(
    payload: dict[str, object], path: Path, parser: argparse.ArgumentParser
) -> None:
    if payload.get("schema") != FIT_SCHEMA:
        parser.error(
            f"{path} uses crash-fit schema {payload.get('schema')!r}; "
            "rerun the cumulative causal fitter instead of loading obsolete parameters"
        )


def _select_fit_cell(
    cells: tuple[CrashFitCell, ...],
    selector: str,
    parser: argparse.ArgumentParser,
) -> CrashFitCell:
    """Resolve one explicit source recording without fuzzy matching."""
    if not cells:
        parser.error("the selected corpus contains no fitting cells")
    if selector.isdecimal():
        index = int(selector) - 1
        if 0 <= index < len(cells):
            return cells[index]
        parser.error(f"--fit-cell index must be between 1 and {len(cells)}")
    matches = [cell for cell in cells if cell.label == selector]
    if len(matches) == 1:
        return matches[0]
    labels = "\n  ".join(cell.label for cell in cells)
    parser.error(f"--fit-cell must exactly match one of:\n  {labels}")


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
    swell_events = tuple(
        CrashEvent(
            0.25 + 0.125 * event_index,
            0.45,
            1.0,
            0.60,
            0x53574580 + event_index,
        )
        for event_index in range(32)
    )
    clips.append(
        AuditionClip(
            "Repeated edge hits, constant strength swell",
            AudioBuffer(render_crash_sequence(fit, 10.0, swell_events), 48000),
        )
    )
    return tuple(clips)


if __name__ == "__main__":
    main()
