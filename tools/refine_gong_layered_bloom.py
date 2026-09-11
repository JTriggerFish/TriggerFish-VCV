"""Keep the user's harmonic core; fit upper bloom with a smooth bar balance.

No centre placement, per-mode damping, extra knots, limiter or normalization.
The low-band target is explicitly the user's envelope minus 16 dB: it remains
about 7 dB above this reference, without dominating the whole instrument.
Upper-band targets remain the reference at its saved gain. This is a hybrid
sound-design objective, not a claim of matching the complete reference.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from scipy.signal import stft

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint


def balanced_parameters(base, low=-16, upper=24):
    """A smooth observation curve; existing low-pair balance stays intact."""
    result = dict(base)
    for i in range(32):
        frequency = base[f"resolved_frequency_{i}"]
        x = np.clip(np.log2(frequency / 300) / np.log2(3000 / 300), 0, 1)
        weight = x * x * (3 - 2 * x)
        level = base[f"resolved_level_{i}"]
        if level > -71.99:
            result[f"resolved_level_{i}"] = float(
                np.clip(level + low * (1 - weight) + upper * weight, -72, 6)
            )
    return result


class LayeredBloomLoss:
    units = "hybrid band-envelope dB error; user core and reference upper bloom"
    bands = [(80, 300), (300, 900), (900, 3000), (3000, 7000), (7000, 14000)]
    times = [0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6]

    def __init__(self, reference, user, rate):
        self.rate, self.frames = rate, len(reference)
        self.target = self.envelopes(reference)
        self.target[0] = self.envelopes(user)[0] - 16
        self.specification = dict(
            version="gong-layered-bloom-v1",
            bands=self.bands,
            times=self.times,
            fft_size=4096,
            hop_seconds=0.01,
            floor_db=-100,
            low_target="user snapshot 80–300 Hz minus 16 dB",
            upper_target="unaltered reference",
            normalization=False,
            weights=dict(low_body=2, middle=1, high=1, high_rise=1),
        )

    def envelopes(self, audio):
        if len(audio) != self.frames or not np.isfinite(audio).all():
            raise ValueError("Expected finite reference-length mono audio")
        f, t, z = stft(
            audio, self.rate, nperseg=4096, noverlap=4096 - round(0.01 * self.rate)
        )
        power = np.array(
            [np.sum(abs(z[(f >= lo) & (f < hi)]) ** 2, axis=0) for lo, hi in self.bands]
        )
        return 10 * np.log10(
            np.maximum(
                1e-10,
                np.array(
                    [
                        power[:, (t >= lo) & (t < hi)].mean(axis=1)
                        for lo, hi in zip(self.times[:-1], self.times[1:])
                    ]
                ).T,
            )
        )

    def diagnostics(self, audio):
        db = self.envelopes(audio)
        error = db - self.target
        rise = error[3:, 3:8] - error[3:, 1:2]
        rms = lambda x: float(np.sqrt(np.mean(x**2)))
        return dict(
            low=rms(error[0]),
            middle=rms(error[1:3]),
            high=rms(error[3:]),
            rise=rms(rise),
            envelopes_db=db.tolist(),
        )

    def score(self, audio):
        d = self.diagnostics(audio)
        return 2 * d["low"] + d["middle"] + d["high"] + d["rise"]


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        fit = json.loads(args.fit.read_text(encoding="utf8"))
        response = r.request(command="renderSnapshot", fit=fit, seconds=6)
        base = {
            k: v
            for node in response["fit"]["instrument"]["nodes"]
            for k, v in node["parameters"].items()
        }
        if response["fit"]["controls"]["event"] != r.metadata["event"]:
            raise ValueError("The source gesture must match the standard gong event")
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "user.fit.json").write_text(
            json.dumps(fit, indent=2), encoding="utf8"
        )
        reference = aligned_reference(r, 6)
        loss = LayeredBloomLoss(reference, r.decode(response["pcm"]), r.sample_rate)
        checkpoint(
            r, loss, args.output / "baseline", "User tuned gong", base, reference, []
        )
        rows = []
        balanced = balanced_parameters(base)
        for exponent in (0.05, 0.2, 0.5, 0.8):
            for sensitivity in (0, 0.3, 0.63):
                for rate in (1, 4, 10, 16):
                    p = dict(
                        balanced,
                        bloom_rate=rate,
                        bloom_energy_acceleration=exponent,
                        bloom_energy_sensitivity=sensitivity,
                    )
                    audio = r.render(p, 6)
                    rows.append(
                        dict(
                            parameters=p,
                            score=loss.score(audio),
                            diagnostics=loss.diagnostics(audio),
                        )
                    )
            print(
                json.dumps(dict(exponent=exponent, best=min(x["score"] for x in rows))),
                flush=True,
            )
        rows.sort(key=lambda x: x["score"])
        (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
        checkpoint(
            r,
            loss,
            args.output / "screen-best",
            "Gong bloom screen",
            rows[0]["parameters"],
            reference,
            [
                dict(
                    stage="coarse dynamics screen",
                    trials=len(rows),
                    specification=loss.specification,
                )
            ],
        )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
