"""Saved-candidate velocity/restrike audit using the actual WASM renderer."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import stft

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def measure(audio, rate):
    if not np.isfinite(audio).all():
        raise ValueError("Nonfinite output")
    f, t, z = stft(audio, rate, nperseg=4096, noverlap=3072)
    power = abs(z) ** 2
    region = (t >= 0.2) & (t < 1.5)
    body = power[(f >= 80) & (f < 800)][:, region].sum()
    high = power[(f >= 5000) & (f < 15000)][:, region].sum()
    return dict(
        sample_peak_db=float(20 * np.log10(max(abs(audio).max(), 1e-20))),
        energy=float(np.sum(audio**2) / rate),
        upper_to_body_db=float(10 * np.log10(max(high, 1e-20) / max(body, 1e-20))),
    )


def run(args):
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(r, args.directory)
        parameters = saved["parameters"]
        output = args.directory / "playability"
        output.mkdir(exist_ok=True)
        report = dict(
            event=r.metadata["event"],
            sample_rate=r.sample_rate,
            velocities=[],
            repeats=[],
        )
        for strength in (0.3, 0.5, 0.76, 1):
            response = r.request(
                command="renderSequence",
                parameters=parameters,
                seconds=6,
                hits=[dict(time=0, strength=strength, seed=1675)],
            )
            audio = r.decode(response["pcm"])
            report["velocities"].append(
                dict(strength=strength, **measure(audio, r.sample_rate))
            )
            write_wav(
                output / f"strength-{strength}.wav", AudioBuffer(audio, r.sample_rate)
            )
        for interval, count in ((0.5, 4), (0.125, 8)):
            hits = [dict(time=i * interval, seed=1675 + i) for i in range(count)]
            audio = r.decode(
                r.request(
                    command="renderSequence",
                    parameters=parameters,
                    seconds=6,
                    hits=hits,
                )["pcm"]
            )
            report["repeats"].append(
                dict(interval=interval, count=count, **measure(audio, r.sample_rate))
            )
            write_wav(
                output / f"repeated-{interval}.wav", AudioBuffer(audio, r.sample_rate)
            )
        energy = [row["energy"] for row in report["velocities"]]
        if not all(b > a for a, b in zip(energy[:-1], energy[1:])):
            raise ValueError("Increasing strength did not increase total output energy")
        (output / "audit.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=2))
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    run(parser.parse_args())
