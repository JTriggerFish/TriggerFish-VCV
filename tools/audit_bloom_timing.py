"""Measure the actual UI timing gesture, without changing any saved fit."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import stft

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def measure(samples, rate):
    frequencies, times, spectrum = stft(
        samples, rate, nperseg=2048, noverlap=1536, boundary="zeros"
    )
    rows = []
    for low, high in ((1000, 3000), (3000, 8000), (8000, 15000)):
        power = np.sum(
            abs(spectrum[(frequencies >= low) & (frequencies < high)]) ** 2, axis=0
        )
        early = times <= 2
        cumulative = np.cumsum(power[early])
        rows.append(
            dict(
                band_hz=[low, high],
                peak_seconds=float(times[np.argmax(power)]),
                median_first_2s=float(
                    times[early][np.searchsorted(cumulative, cumulative[-1] / 2)]
                ),
                initial_fraction=float(
                    power[times < 0.1].sum() / max(power.sum(), 1e-30)
                ),
            )
        )
    return rows


def run(args):
    if args.source is not None and len(args.targets) != 1:
        raise ValueError("Use one target when auditing a saved checkpoint")
    result = {}
    for name in args.targets:
        renderer = WorkbenchRenderer(
            os.environ["EMSDK_NODE"], name + "-standard", Path.cwd()
        )
        try:
            baseline = (
                renderer.initial
                if args.source is None
                else verify_candidate(renderer, args.source)["parameters"]
            )
            rows = []
            for position in (-1, -0.5, 0, 0.5, 1):
                expansion = renderer.request(
                    command="bloomTiming",
                    parameters=baseline,
                    position=position,
                )
                parameters = dict(baseline, **expansion["values"])
                audio = renderer.render(parameters, 6)
                if not np.isfinite(audio).all():
                    raise ValueError("Nonfinite timing render")
                rows.append(
                    dict(
                        position=position,
                        **expansion,
                        measurements=measure(audio, renderer.sample_rate),
                    )
                )
            result[name] = dict(
                parameters=baseline,
                renderer=renderer.metadata,
                reference=measure(aligned_reference(renderer, 6), renderer.sample_rate),
                rows=rows,
            )
        finally:
            renderer.close()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", nargs="+", default=["crash", "gong"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--source",
        type=Path,
        help="Audit a validated checkpoint instead of the published preset",
    )
    run(parser.parse_args())
