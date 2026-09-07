"""Independent before/reference/candidate checks for the reported kick defects."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from triggerfish_percussion.transforms import StftConfig, stft
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.fit_rerender import load_checked_start

REGIONS = ((0, 0.08), (0.08, 0.16), (0.16, 0.26), (0.26, 0.4))
BANDS = ((20, 45), (48, 65), (78, 100), (100, 200), (540, 650), (760, 860))


def measurements(audio, rate):
    value = stft(audio, rate, StftConfig(8192, 128))
    rows = []
    for lo, hi in BANDS:
        bins = (value.frequencies_hz >= lo) & (value.frequencies_hz < hi)
        rows.append(
            [
                value.power[bins][
                    :, (value.times_seconds >= a) & (value.times_seconds < b)
                ]
                .sum(axis=0)
                .mean()
                for a, b in REGIONS
            ]
        )
    return np.maximum(rows, 1e-20)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        required=True,
        help="Frozen baseline search.json; never the live preset",
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        action="append",
        required=True,
        help="Explicit candidate search.json (repeatable)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New audit JSON; existing reports are never overwritten",
    )
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Audit already exists; choose a new output path")
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", root)
    try:
        candidates, sources = [], []
        for index, path in enumerate([args.baseline, *args.candidate]):
            saved, provenance = load_checked_start(renderer, path)
            name = "baseline" if index == 0 else f"candidate-{index}"
            candidates.append((name, saved["parameters"]))
            sources.append(dict(name=name, **provenance))
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset:]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        target = measurements(reference, rate)
        rows = []
        for name, values in candidates:
            for seed in (1449, 1451, 1452):
                audio = renderer.render(values, 1.2, seed)
                power = measurements(audio, rate)
                error = 10 * np.log10(power / target)
                rows.append(
                    dict(
                        name=name,
                        seed=seed,
                        band_error_db=error.tolist(),
                        peak=float(np.max(np.abs(audio))),
                        decay_ratio_error_db=(error[:4, 3] - error[:4, 0]).tolist(),
                    )
                )
        report = dict(
            bands_hz=BANDS,
            regions_seconds=REGIONS,
            results=rows,
            sources=sources,
            renderer=renderer.metadata,
            reference_db=(10 * np.log10(target)).tolist(),
            reference_peak=float(np.max(np.abs(reference))),
        )
        with args.output.open("x", encoding="utf8") as stream:
            json.dump(report, stream, indent=2)
        for row in rows:
            print(json.dumps(row))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
