"""Small audible-rate/clarity screen; no independent upper-mode fitting."""

import json
import os
import argparse
from pathlib import Path
import numpy as np
import torch
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.low_mode_beating import LowModeBeating
from triggerfish_percussion.workbench_fit_baseline import check_reference
from crash_beating_common import CrashObjective
from fit_stretched_gong import checkpoint


def main(source_path):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    output = Path("build/crash-paired-ring")
    try:
        source = json.loads(source_path.read_text(encoding="utf8"))
        check_reference(source["metadata"], renderer.metadata)
        base = source["parameters"]
        if set(base) != set(renderer.initial) or any(
            not np.isfinite(base[d["key"]])
            or not d["minimum"] <= base[d["key"]] <= d["maximum"]
            for d in renderer.metadata["descriptors"]
        ):
            raise ValueError(
                "Starting parameters must match the current public surface"
            )
        reference = aligned_reference(renderer, 6)
        objective = CrashObjective(reference, renderer.sample_rate, 60)
        beating = LowModeBeating(reference, renderer.sample_rate)
        rows = []
        for rate in (0.75, 1.25, 1.5, 2):
            for clarity in (0, 0.4, 0.7, 1):
                p = dict(
                    base,
                    field_distribution=3,
                    field_doublet_split=rate,
                    resolved_frequency_0=125.9,
                )
                # One smooth low-end clarity edit, baked into existing handles.
                for i in range(32):
                    f = p[f"resolved_frequency_{i}"]
                    weight = np.clip(np.log(900 / f) / np.log(900 / 200), 0, 1)
                    p[f"resolved_turbulence_{i}"] *= 1 + (clarity - 1) * weight
                audio = renderer.render(p, 6)
                bands = beating.analyze(audio)
                c = objective.components(audio)
                low = beating.score_rows(bands)
                row = dict(
                    rate=rate,
                    clarity=clarity,
                    parameters=p,
                    components=c,
                    low_score=low,
                    low_band=bands[0],
                    score=objective.score_components(c) + 0.3 * low,
                )
                rows.append(row)
                print(
                    json.dumps(
                        dict(
                            rate=rate,
                            clarity=clarity,
                            score=row["score"],
                            low_score=low,
                            pulse=bands[0]["dominant_hz"],
                            fast=bands[0]["fast_fraction"],
                            power=bands[0]["power"],
                        )
                    ),
                    flush=True,
                )
        (output / "screen.json").write_text(json.dumps(rows, indent=2))
        # The coarse modulation bands alone cannot distinguish 1.25 from 2 Hz.
        # Constrain this first trial by the independently measured slow pulse,
        # then compare low-band texture. Do not call the composite minimum a
        # better fit when it visibly misses the user's requested mechanism.
        pulse = beating.target[0]["dominant_hz"]
        nearest = min(abs(r["rate"] - pulse) for r in rows)
        eligible = [r for r in rows if abs(r["rate"] - pulse) <= nearest + 1e-6]
        best = min(eligible, key=lambda r: r["low_score"])
        checkpoint(
            renderer,
            objective,
            output / "candidate",
            "Crash — paired ring trial",
            best["parameters"],
            reference,
            [
                dict(
                    stage="paired ring screen",
                    source_parameters=str(source_path),
                    candidates=[
                        {
                            k: v
                            for k, v in r.items()
                            if k not in ("parameters", "low_band")
                        }
                        for r in rows
                    ],
                    selection="nearest measured low pulse, then low-band texture; all other metrics reported separately",
                )
            ],
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("build/crash-paired-ring/before/search.json"),
    )
    main(parser.parse_args().source)
