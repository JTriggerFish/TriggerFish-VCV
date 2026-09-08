"""Small crash velocity-colour screen with a fixed nominal-strike guard.

Fit 48/72/96; report but do not select on 24/120. Strength stays linear and each
reference keeps the same source-family gain. This fits one continuous patch,
not independent per-velocity gain corrections or an interpolated preset grid.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def run(args):
    torch.set_num_threads(1)
    renderers = {}
    try:
        for velocity in (24, 48, 72, 96, 120):
            renderers[velocity] = WorkbenchRenderer(
                os.environ["EMSDK_NODE"],
                "crash-standard",
                Path.cwd(),
                cell=dict(velocity=velocity),
            )
        standard = renderers[72]
        saved = verify_candidate(standard, args.source)
        seconds = saved["duration_seconds"]
        references = {v: aligned_reference(r, seconds) for v, r in renderers.items()}
        losses = {
            v: AuralossMel(references[v], r.sample_rate) for v, r in renderers.items()
        }
        shape = SpectralBloomLoss(references[72], standard.sample_rate)
        initial = saved["parameters"]
        candidates = [("baseline", initial)]
        for colour in (0, 1, 2, 3, 6, 8):
            # First-order compensation around the fixed nominal strike, using
            # the exposed initial tilt. No runtime macro or hidden coefficient.
            tilt = initial["body_brightness"] + (
                initial["velocity_brightness"] - colour
            ) * (standard.metadata["event"]["strength"] - 0.8)
            candidates.append(
                (
                    f"velocity-colour-{colour}",
                    dict(
                        initial,
                        velocity_brightness=colour,
                        body_brightness=float(np.clip(tilt, -72, 24)),
                    ),
                )
            )
        rows, best, chosen = [], float("inf"), initial
        nominal_limit = shape_limit = None
        for name, parameters in candidates:
            audios = {v: r.render(parameters, seconds) for v, r in renderers.items()}
            scores = {v: losses[v].score(x) for v, x in audios.items()}
            bloom = float(np.linalg.norm(shape.residual(audios[72])))
            if nominal_limit is None:
                nominal_limit, shape_limit = scores[72] + 0.05, bloom + 0.25
            eligible = scores[72] <= nominal_limit and bloom <= shape_limit
            score = (scores[48] + 2 * scores[72] + scores[96]) / 4
            if eligible and score < best:
                best, chosen = score, parameters
            row = dict(
                name=name,
                scores=scores,
                bloom=bloom,
                eligible=eligible,
                training_score=score,
                parameters=parameters,
            )
            rows.append(row)
            print(
                json.dumps({k: v for k, v in row.items() if k != "parameters"}),
                flush=True,
            )
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            standard,
            losses[72],
            args.output,
            seconds,
            "Crash - bloom and velocity colour",
        )
        search.parameters = chosen
        search.history.append(
            dict(
                parent=str(args.source.resolve()),
                stage="velocity-colour screen",
                training_velocities=[48, 72, 96],
                holdout_velocities=[24, 120],
                training_weights=[1, 2, 1],
                nominal_mel_limit=nominal_limit,
                nominal_bloom_limit=shape_limit,
                candidates=rows,
                reference_metadata={v: r.metadata for v, r in renderers.items()},
            )
        )
        search.save()
        write_wav(
            args.output / "reference.wav",
            AudioBuffer(references[72], standard.sample_rate),
        )
        verify_candidate(standard, args.output)
    finally:
        for renderer in renderers.values():
            renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
