"""Fit only the root and stretch of an existing harmonic grid, in actual Wasm.

The bounded coordinates are root Hz and highest-mode Hz, avoiding invalid
combinations above 15 kHz. Both are translated by the shared UI generator.
The objective includes initial 400-ms Mel, full Mel and a bloom guard.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from refine_gong_strike import StrikeLoss
from triggerfish_percussion.attack_ridge_loss import AttackRidgeLoss


class FineStrikeLoss(StrikeLoss):
    def __init__(self, reference, rate, baseline):
        super().__init__(reference, rate, baseline)
        self.ridges = AttackRidgeLoss(reference, rate)
        self.specification.update(
            attack_ridges=self.ridges.specification,
            weights=dict(early_mel=0.5, early_linear_stft=0.5, full_mel=1),
        )

    def score(self, samples):
        excess = max(0.0, np.linalg.norm(self.bloom.residual(samples)) - self.limit)
        return (
            0.5 * self.early.score(samples[: self.frames])
            + 0.5 * self.ridges.score(samples)
            + self.full.score(samples)
            + excess**2
        )


def run(args):
    torch.set_num_threads(1)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        saved = verify_candidate(renderer, args.source)
        reference = aligned_reference(renderer, 6)
        baseline = renderer.render(saved["parameters"], 6)
        loss_type = FineStrikeLoss if args.fine_ridges else StrikeLoss
        loss = loss_type(reference, renderer.sample_rate, baseline)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(
            renderer, loss, args.output, 6, "Gong - stretched pitch refinement"
        )
        original = saved["parameters"]
        count = sum(original[f"resolved_level_{i}"] > -71.99 for i in range(32))
        root = original["resolved_frequency_0"]
        top = original[f"resolved_frequency_{count-1}"]
        initial_stretch = renderer.request(
            command="modalTemplateStretch",
            settings=dict(
                fundamental=root,
                count=count,
                topFrequency=top,
                harmonicCore=args.harmonic_core,
            ),
        )["stretch"]
        initial_points = renderer.request(
            command="modalTemplate",
            settings=dict(
                family="harmonic",
                fundamental=root,
                count=count,
                stretch=initial_stretch,
                harmonicCore=args.harmonic_core,
                minimumFrequency=1,
            ),
        )["points"]
        if not all(
            np.isclose(
                original[f"resolved_frequency_{i}"],
                initial_points[i]["frequency"],
            )
            for i in range(count)
        ):
            raise ValueError(
                "Source must use the current protected-core stretch law; historical power-law fits cannot be resumed here"
            )
        low, high = np.array([75.0, 7500.0]), np.array([145.0, 15000.0])
        best = [loss.score(baseline), dict(original), [root, top]]
        trials = []

        def objective(unit):
            root_hz, top_hz = low + np.asarray(unit) * (high - low)
            try:
                stretch = renderer.request(
                    command="modalTemplateStretch",
                    settings=dict(
                        fundamental=float(root_hz),
                        count=count,
                        topFrequency=float(top_hz),
                        harmonicCore=args.harmonic_core,
                    ),
                )["stretch"]
            except RuntimeError as error:
                if "Top frequency is outside" not in str(error):
                    raise
                trials.append(
                    dict(root=float(root_hz), top=float(top_hz), rejected=str(error))
                )
                return 1e6
            points = renderer.request(
                command="modalTemplate",
                settings=dict(
                    family="harmonic",
                    fundamental=float(root_hz),
                    count=count,
                    stretch=stretch,
                    harmonicCore=args.harmonic_core,
                    minimumFrequency=1,
                    maximumFrequency=15000,
                ),
            )["points"]
            parameters = dict(
                original,
                **{
                    f"resolved_frequency_{i}": p["frequency"]
                    for i, p in enumerate(points)
                },
            )
            score = loss.score(search.audio(parameters))
            trials.append(
                dict(
                    root=float(root_hz), top=float(top_hz), stretch=stretch, score=score
                )
            )
            if score < best[0]:
                best[:] = [score, parameters, [float(root_hz), float(top_hz)]]
            if len(trials) % 10 == 0:
                print(
                    json.dumps(
                        dict(evaluations=len(trials), best=best[0], geometry=best[2])
                    ),
                    flush=True,
                )
            return score

        minimize(
            objective,
            (np.array([root, top]) - low) / (high - low),
            method="Powell",
            bounds=[(0, 1), (0, 1)],
            options=dict(maxfev=140, xtol=0.001, ftol=0.0002),
        )
        search.parameters = best[1]
        search.history.append(
            dict(
                parent=str(args.source.resolve()),
                stage="root and stretch",
                method="bounded Powell; root/top coordinates",
                bounds=[low.tolist(), high.tolist()],
                before=loss.score(baseline),
                after=best[0],
                trials=trials,
                fixed_parameters={
                    k: v
                    for k, v in original.items()
                    if not k.startswith("resolved_frequency_")
                },
            )
        )
        search.save()
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fine-ridges", action="store_true")
    parser.add_argument("--harmonic-core", type=int, choices=range(1, 9), default=4)
    run(parser.parse_args())
