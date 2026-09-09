"""Screen texture-only gong edits; preserve its modal series, bloom and damping.

Use reference modulation and absolute spectral error together. A lower
modulation score alone is not a reason to replace a useful instrument fit.
"""

import json
import os
from pathlib import Path

import numpy as np
import torch
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.low_mode_beating import LowModeBeating
from triggerfish_percussion.modal_texture_loss import ModalTextureLoss
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def variants(base):
    yield "before", base
    for rate in (0.75, 1.5, 3):
        yield f"doublets-{rate}", dict(base, field_doublet_split=rate)
    for blur in (0.012, 0.025):
        yield f"doublets-3-blur-{blur}", dict(
            base, field_doublet_split=3, field_phase_bandwidth=blur
        )
    for slope in (0.3, 0.4, 0.5):
        yield f"clean-low-{slope}", dict(
            base, field_doublet_split=3, field_turbulence_slope=slope
        )
    for rate in (0.75, 1.25):
        for depth in (0, 0.15, 0.3, 0.45):
            yield f"paired-{rate}-{depth}", dict(
                base,
                field_distribution=3,
                field_doublet_split=rate,
                field_beat_depth=depth,
                field_beat_rate_tilt=0.25,
            )


def main():
    torch.set_num_threads(1)
    output = Path("build/gong-beating-recovery")
    output.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        reference = aligned_reference(renderer, 6)
        beating = LowModeBeating(reference, renderer.sample_rate)
        mel = ReferenceFloorMel(reference, renderer.sample_rate, 60)
        texture = ModalTextureLoss(reference, renderer.sample_rate)
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        records = []
        seed = renderer.metadata["event"]["seed"]
        for name, parameters in variants(renderer.initial):
            rows = []
            for offset in (0, 911, 1601):
                audio = renderer.render(parameters, 6, seed + offset)
                bands = beating.analyze(audio)
                rows.append(
                    dict(
                        seed=seed + offset,
                        mel=mel.score(audio),
                        texture=texture.score(audio),
                        beating=beating.score_rows(bands),
                        power=[b["power"] for b in bands],
                    )
                )
                if offset == 0:
                    write_wav(
                        output / f"{name}.wav", AudioBuffer(audio, renderer.sample_rate)
                    )
            row = dict(
                name=name,
                parameters=parameters,
                seeds=rows,
                mel=float(np.mean([r["mel"] for r in rows])),
                beating=float(np.mean([r["beating"] for r in rows])),
            )
            records.append(row)
            print(
                json.dumps({k: row[k] for k in ("name", "mel", "beating")}), flush=True
            )
        (output / "screen.json").write_text(
            json.dumps(
                dict(
                    reference=beating.target,
                    metadata=renderer.metadata,
                    variants=records,
                ),
                indent=2,
            ),
            encoding="utf8",
        )
    finally:
        renderer.close()


if __name__ == "__main__":
    main()
