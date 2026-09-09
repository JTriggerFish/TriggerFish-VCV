"""Separate low-mode spacing, stochastic phase blur and energy diffusion."""

import json
import os
from pathlib import Path
import numpy as np
import torch
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.low_mode_beating import LowModeBeating
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.audio_io import AudioBuffer, write_wav


def variants(base):
    yield "current", dict(base)
    yield "no-blur", dict(base, field_phase_bandwidth=0)
    yield "no-diffusion", dict(base, bloom_rate=0)
    yield "neither", dict(base, field_phase_bandwidth=0, bloom_rate=0)
    for scale in (0.25, 0.5, 0.75):
        p = dict(base)
        for i in range(32):
            f = p[f"resolved_frequency_{i}"]
            weight = np.clip(np.log(1000 / f) / np.log(1000 / 450), 0, 1)
            p[f"resolved_turbulence_{i}"] *= 1 + (scale - 1) * weight
        yield f"narrow-lows-{scale}", p
    p = dict(base)
    for i in range(32):
        if p[f"resolved_frequency_{i}"] < 600:
            p[f"resolved_allocation_{i}"] = 3
    yield "denser-lows", p
    for split in (0.75, 1.5, 3):
        yield f"doublets-{split}", dict(
            base, field_distribution=2, field_doublet_split=split
        )
    # Test a coherent low pair with ordinary visible handles before changing
    # the engine. Only this low pair is placed individually; the upper series
    # and its texture stay intact. Never overwrite an active user handle.
    free = next((i for i in range(32) if base[f"resolved_level_{i}"] <= -71.99), None)
    if free is not None:
        lowest = min(
            (i for i in range(32) if base[f"resolved_level_{i}"] > -71.99),
            key=lambda i: base[f"resolved_frequency_{i}"],
        )
        for split in (0.75, 1.25, 2):
            p = dict(base)
            p[f"resolved_frequency_{free}"] = p[f"resolved_frequency_{lowest}"] + split
            p[f"resolved_level_{free}"] = p[f"resolved_level_{lowest}"] - 3
            p[f"resolved_turbulence_{lowest}"] = 0
            p[f"resolved_turbulence_{free}"] = 0
            yield f"paired-low-{split}", p
        # Five-second reference spectrum resolves prominent peaks near these
        # frequencies (0.2-Hz bins). An explicit target-specific experiment,
        # not a new hard-coded engine constant or general cymbal tuning law.
        for balance in (-6, -8, -10):
            p = dict(base)
            p[f"resolved_frequency_{lowest}"] = 126.6
            p[f"resolved_frequency_{free}"] = 125.2
            p[f"resolved_level_{free}"] = p[f"resolved_level_{lowest}"] + balance
            p[f"resolved_turbulence_{lowest}"] = 0
            p[f"resolved_turbulence_{free}"] = 0
            yield f"measured-pair-{abs(balance)}", p


def main():
    torch.set_num_threads(1)
    output = Path("build/crash-low-beating")
    output.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "crash-standard", Path.cwd())
    try:
        reference = aligned_reference(renderer, 6)
        loss = LowModeBeating(reference, renderer.sample_rate)
        spectral = ReferenceFloorMel(reference, renderer.sample_rate, 60)
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        results = dict(
            reference=loss.target, specification=loss.specification, variants=[]
        )
        for name, p in variants(renderer.initial):
            audio = renderer.render(p, 6)
            rows = loss.analyze(audio)
            row = dict(
                name=name,
                low_score=loss.score_rows(rows),
                spectral_score=spectral.score(audio),
                bands=rows,
                parameters=p,
            )
            results["variants"].append(row)
            write_wav(output / f"{name}.wav", AudioBuffer(audio, renderer.sample_rate))
            print(
                json.dumps(
                    dict(
                        name=name,
                        low_score=row["low_score"],
                        spectral_score=row["spectral_score"],
                        fast=[r["fast_fraction"] for r in rows],
                        dominant=[r["dominant_hz"] for r in rows],
                    )
                ),
                flush=True,
            )
        (output / "audit.json").write_text(json.dumps(results, indent=2))
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
