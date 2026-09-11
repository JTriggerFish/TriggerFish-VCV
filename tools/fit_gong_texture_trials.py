"""EQ-free texture-family screens through the real WASM, not a new synthesizer."""

import argparse
import json
import os
from itertools import product
from pathlib import Path

import numpy as np

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from gong_texture_comparison import GongTextureComparison


def cloud(base, centre, width):
    p = dict(
        base, field_distribution=4, field_satellite_density=1, field_packet_spread=width
    )
    high = [i for i in range(32) if p[f"resolved_frequency_{i}"] >= 3000]
    for i in high:
        p[f"resolved_level_{i}"] = -72
    anchor = high[0]
    p[f"resolved_frequency_{anchor}"] = centre
    p[f"resolved_level_{anchor}"] = -12
    p[f"resolved_turbulence_{anchor}"] = 0.65
    p[f"resolved_allocation_{anchor}"] = 4
    return p


def trials(base):
    for blur, tilt, spread in product((0.025, 0.075, 0.2), (0, 0.5, 1), (1.8, 3.5)):
        yield "blur", dict(
            base,
            field_motion_depth=0,
            field_phase_bandwidth=blur,
            field_phase_tilt=tilt,
            field_packet_spread=spread,
            field_satellite_density=1,
        )
    for depth, speed, spread in product((0.6, 1.5, 3), (40, 120, 200), (1.8, 3.5)):
        yield "movement", dict(
            base,
            field_motion_depth=depth,
            field_motion_rate=speed,
            field_motion_sharing=0.15,
            field_phase_bandwidth=0,
            field_packet_spread=spread,
            field_satellite_density=1,
        )
    for centre, width, mechanism, rate in product(
        (7500, 9500), (2.5, 4), ("blur", "movement"), (4, 10)
    ):
        p = cloud(base, centre, width)
        p.update(
            bloom_rate=rate,
            field_phase_bandwidth=0.05 if mechanism == "blur" else 0,
            field_phase_tilt=0.5,
            field_motion_depth=1.5 if mechanism == "movement" else 0,
            field_motion_rate=120,
            field_motion_sharing=0.15,
        )
        yield "cloud", p


def run(args):
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        base = dict(r.initial, output_eq_enabled=0)
        body = r.render(base, 6, 1675)
        ref = aligned_reference(r, 6)
        loss = GongTextureComparison(ref, body, r.sample_rate)
        args.output.mkdir(parents=True, exist_ok=True)
        initial = dict(
            parameters=base, metadata=r.metadata, specification=loss.specification
        )
        (args.output / "baseline.json").write_text(json.dumps(initial, indent=2))
        base["field_wander_hz"] = 0  # Keep the two movement mechanisms separable.
        rows = []
        for index, (family, p) in enumerate(trials(base)):
            m = loss.measure(r.render(p, 6, 1675))
            # Screen dynamics/texture; amplitude polishing follows separately.
            score = (
                m["upper_shape_db"]
                + 0.3 * m["upper_db"]
                + 0.4 * m["body_db"]
                + m["ridge_error_db"]
            )
            rows.append(dict(family=family, parameters=p, metrics=m, score=score))
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            if index % 6 == 0:
                print(
                    json.dumps(dict(index=index, family=family, score=score)),
                    flush=True,
                )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
