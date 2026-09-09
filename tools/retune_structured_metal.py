"""Move the series root while preserving its top and broad prominence curve."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from fit_structured_metal_texture import Objective


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], args.target + "-standard", Path.cwd()
    )
    try:
        source = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        objective = Objective(reference, renderer.sample_rate)
        original = source["parameters"]
        active = [i for i in range(32) if original[f"resolved_level_{i}"] > -71.99]
        frequencies = np.array([original[f"resolved_frequency_{i}"] for i in active])
        amplitudes = 10 ** (
            np.array([original[f"resolved_level_{i}"] for i in active]) / 20
        )
        knots = (120, 600, 3000, 15000)
        weights = interpolation_weights(frequencies, knots)
        curve = np.linalg.lstsq(weights, amplitudes, rcond=None)[0]
        if np.linalg.norm(weights @ curve - amplitudes) > 1e-6 * np.linalg.norm(
            amplitudes
        ):
            raise ValueError("Source must have four-coordinate prominence")
        seeds = [renderer.metadata["event"]["seed"] + offset for offset in (0, 101)]

        def score(p):
            return float(
                np.mean([objective.score(renderer.render(p, 6, s)) for s in seeds])
            )

        best, best_score = original, score(original)
        rows = [dict(root=float(frequencies[0]), score=best_score, baseline=True)]
        for root in args.roots:
            settings = dict(
                fundamental=root,
                count=len(active),
                harmonicCore=4,
                topFrequency=float(frequencies[-1]),
            )
            stretch = renderer.request(
                command="modalTemplateStretch", settings=settings
            )["stretch"]
            points = renderer.request(
                command="modalTemplate",
                settings=dict(
                    family="harmonic",
                    fundamental=root,
                    count=len(active),
                    harmonicCore=4,
                    stretch=stretch,
                    minimumFrequency=1,
                ),
            )["points"]
            f = np.array([point["frequency"] for point in points])
            levels = 20 * np.log10(interpolation_weights(f, knots) @ curve)
            p = dict(original)
            for i, frequency, level in zip(active, f, levels):
                p[f"resolved_frequency_{i}"] = float(frequency)
                p[f"resolved_level_{i}"] = float(np.clip(level, -72, 6))
            for centre in args.allocation_centres:
                for split in args.splits:
                    trial = dict(p, field_doublet_split=split)
                    # A broad log-frequency allocation ramp, not independent
                    # local fits: zero below half-centre, one above twice-centre.
                    if centre > 0:
                        for i, frequency in zip(active, f):
                            x = float(
                                np.clip(0.5 + 0.5 * np.log2(frequency / centre), 0, 1)
                            )
                            trial[f"resolved_allocation_{i}"] = x * x * (3 - 2 * x)
                    value = score(trial)
                    row = dict(
                        root=root,
                        stretch=stretch,
                        allocation_centre=centre,
                        split=split,
                        score=value,
                    )
                    rows.append(row)
                    print(json.dumps(row), flush=True)
                    if value < best_score:
                        best, best_score = trial, value
        search = checkpoint(
            renderer,
            objective,
            args.output,
            args.target.title() + " — structured beating fit",
            best,
            reference,
            source["history"]
            + [dict(stage="root screen, top fixed", rows=rows, seeds=seeds)],
        )
        verify_candidate(renderer, search.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["gong", "crash"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--roots", type=float, nargs="+", default=[112, 116, 120, 124, 128, 132]
    )
    parser.add_argument("--allocation-centres", type=float, nargs="+", default=[0])
    parser.add_argument("--splits", type=float, nargs="+", default=[6])
    run(parser.parse_args())
