"""Check the timing gesture in the current WASM, without refitting presets.

Compare high-band energy timing at earlier/centre/later settings. These are
diagnostics, not a universal monotonicity claim for arbitrary sound designs.
Hold decay is deliberately absent so only the three macro coordinates move.
"""

import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import stft

from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def high_band_timing(audio, rate):
    frequency, time, spectrum = stft(audio, rate, nperseg=2048, noverlap=1536)
    selected = (time >= 0) & (time <= 3)
    power = np.sum(
        abs(spectrum[(frequency >= 3000) & (frequency < 12000)]) ** 2, axis=0
    )[selected]
    time = time[selected]
    if not np.isfinite(power).all() or power.sum() <= 0:
        raise ValueError("Invalid high-band energy")
    cumulative = np.cumsum(power) / power.sum()
    return dict(
        peak_seconds=float(time[np.argmax(power)]),
        half_energy_seconds=float(np.interp(0.5, cumulative, time)),
        centroid_seconds=float(np.sum(time * power) / power.sum()),
    )


def main():
    rows = []
    for target in ("gong-standard", "crash-standard"):
        renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], target, Path.cwd())
        try:
            for sensitivity in (renderer.initial["bloom_energy_sensitivity"], 1):
                base = dict(renderer.initial, bloom_energy_sensitivity=sensitivity)
                for position in (-1, 0, 1):
                    changes = renderer.request(
                        command="bloomTiming", parameters=base, position=position
                    )
                    parameters = dict(base, **changes["values"])
                    audio = renderer.render(parameters, 6)
                    row = dict(
                        target=target,
                        sensitivity=sensitivity,
                        position=position,
                        changes=changes,
                        **high_band_timing(audio, renderer.sample_rate),
                    )
                    rows.append(row)
                    print(json.dumps(row), flush=True)
        finally:
            renderer.close()
    directory = Path("build/bloom-timing-ui")
    directory.mkdir(exist_ok=True)
    (directory / "timing-check.json").write_text(
        json.dumps(rows, indent=2), encoding="utf8"
    )


if __name__ == "__main__":
    main()
