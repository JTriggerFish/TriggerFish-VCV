"""Retune shared damping/transport after changing the texture mechanism."""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from gong_texture_comparison import GongTextureComparison
from polish_gong_texture_trials import polish


def variations(p, family, timing_only=False):
    yield dict(p)
    if timing_only:
        for rate, concentration in product((0.75, 1.5, 2.5), (0.02, 0.05, 0.12)):
            yield dict(p, bloom_rate=rate, bloom_energy_acceleration=concentration)
        return
    if family == "cloud":
        anchor = next(
            i
            for i in range(32)
            if p[f"resolved_frequency_{i}"] >= 3000
            and p[f"resolved_level_{i}"] > -71.99
        )
        for centre, width, rate in product((8000, 9000), (5, 7), (4, 10, 16)):
            trial = dict(p, field_packet_spread=width, bloom_rate=rate)
            trial[f"resolved_frequency_{anchor}"] = centre
            for i in range(32):
                trial[f"resolved_allocation_{i}"] = 4 if i == anchor else 0.2
            yield trial
    else:
        for rate, decay in product((2.5, 4, 6), (1.4, 2.14, 3)):
            yield dict(p, bloom_rate=rate, body_decay_seconds_7=decay)


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = json.loads((args.directory / "baseline.json").read_text())["parameters"]
        target = GongTextureComparison(
            aligned_reference(r, 6), r.render(base, 6, 1675), r.sample_rate
        )
        for family in args.families:
            source = f"{family}-final" if args.timing_only else family
            p = json.loads((args.directory / source / "search.json").read_text())[
                "parameters"
            ]
            rows = []
            choices = list(variations(p, family, args.timing_only))
            descriptors = {d["key"]: d for d in r.metadata["descriptors"]}
            for trial in choices:
                for key, value in trial.items():
                    d = descriptors[key]
                    if not d["minimum"] <= value <= d["maximum"]:
                        raise ValueError(
                            f"Trial outside the visible control range: {key}={value}"
                        )
            for trial in choices:
                m = target.measure(r.render(trial, 6, 1675))
                score = (
                    m["upper_shape_db"]
                    + 0.4 * m["upper_db"]
                    + 0.4 * m["body_db"]
                    + 0.5 * m["ridge_error_db"]
                )
                rows.append(dict(parameters=trial, metrics=m, score=score))
                print(
                    json.dumps(dict(family=family, trial=len(rows), score=score)),
                    flush=True,
                )
            (args.directory / f"{family}-dynamics.json").write_text(
                json.dumps(rows, indent=2)
            )
            best = min(rows, key=lambda row: row["score"])
            polish(r, target, best["parameters"], args.directory / f"{family}-final")
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--timing-only", action="store_true")
    parser.add_argument("--families", nargs="+", default=["movement"])
    run(parser.parse_args())
