"""Controlled faults and matched fits. Never publishes a workbench preset."""

import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import time

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.band_region_audit import BandRegionAudit
from triggerfish_percussion.region_fit_loss import RegionFitLoss
from triggerfish_percussion.perceptual_fit_losses import AuralossMel, JtfsLoss
from triggerfish_percussion.remote_perceptual_loss import RemoteJtfsLoss
from triggerfish_percussion.loss_perturbations import fault_examples
from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.kick_quality_checks import check_kick_candidate
from kick_loss_experiment_setup import experiment_starts

ROOT = Path(__file__).resolve().parents[1]


class RegionScalar(RegionFitLoss):
    def score(self, samples):
        return float(np.sum(self.residual(samples) ** 2))


def create_losses(reference, rate, remote, stack):
    import torch

    torch.set_num_threads(2)
    if remote:
        jtfs = RemoteJtfsLoss(
            reference,
            rate,
            "MLBox",
            "/tmp/tf-kick-perceptual-20260905",
            "/home/jt/TriggerFish-VCV/.venv/bin/python",
        )
        stack.callback(jtfs.close)
    else:
        jtfs = JtfsLoss(reference, rate)
    return dict(
        region=RegionScalar(BandRegionAudit(reference, rate), reference, rate),
        mel=AuralossMel(reference, rate),
        mel_a=AuralossMel(reference, rate, True),
        jtfs=jtfs,
    )


def audit_losses(losses, reference, rate, directory):
    rows = []
    examples = [("identity", 0, reference), *fault_examples(reference, rate)]
    for name, severity, audio in examples:
        scores, timings = {}, {}
        for key, loss in losses.items():
            start = time.perf_counter()
            scores[key] = loss.score(audio)
            timings[key] = time.perf_counter() - start
        rows.append(dict(fault=name, severity=severity, scores=scores, seconds=timings))
        print(json.dumps(rows[-1]), flush=True)
    report = dict(rows=rows, objectives={k: v.specification for k, v in losses.items()})
    (directory / "fault-audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf8"
    )
    return report


def score_baselines(renderer, losses, starts, output):
    """Compare to the published patch too, not only deliberately changed starts."""
    rows = []
    for name, parameters in dict(published=renderer.initial, **starts).items():
        audio = renderer.render(parameters, 1.2)
        write_wav(output / f"{name}.wav", AudioBuffer(audio, renderer.sample_rate))
        rows.append(
            dict(
                name=name,
                parameters=parameters,
                scores={k: v.score(audio) for k, v in losses.items()},
            )
        )
    (output / "baselines.json").write_text(json.dumps(rows, indent=2), encoding="utf8")


def run_trials(renderer, reference, losses, starts, bounds, output, budget):
    seed = renderer.metadata["event"]["seed"]
    rows = []
    for name, loss in losses.items():
        for start_name, initial in starts.items():
            label = f"{name}-{start_name}"
            directory = output / label
            directory.mkdir(exist_ok=True)
            search = Search(renderer, loss, directory, 1.2, label, (seed, seed + 11))
            search.parameters = dict(initial)
            write_wav(
                directory / "reference.wav",
                AudioBuffer(reference, renderer.sample_rate),
            )
            record = refine_scalar(search, bounds, budget)
            audio = search.audio(search.parameters, seed)
            cross_scores = {key: value.score(audio) for key, value in losses.items()}
            checks = check_kick_candidate(directory, renderer=renderer)
            rows.append(
                dict(name=label, fit=record, cross_scores=cross_scores, checks=checks)
            )
            (output / "trials.json").write_text(
                json.dumps(rows, indent=2), encoding="utf8"
            )
            print(
                json.dumps(
                    dict(
                        completed=label,
                        before=record["before"],
                        after=record["after"],
                        eligible=checks["eligible"],
                    )
                ),
                flush=True,
            )


def save_manifest(directory, manifest):
    """An audit rerun must never relabel fits from a different experiment."""
    path = directory / "experiment.json"
    payload = json.dumps(manifest, indent=2)
    if path.exists() and json.loads(path.read_text(encoding="utf8")) != json.loads(
        payload
    ):
        raise ValueError("Experiment differs: choose a new output directory")
    path.write_text(payload, encoding="utf8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fit", action="store_true", default=os.environ.get("TF_KICK_LOSS_FIT") == "1"
    )
    parser.add_argument(
        "--remote-jtfs",
        action="store_true",
        default=os.environ.get("TF_KICK_LOSS_REMOTE") == "1",
    )
    parser.add_argument(
        "--budget", type=int, default=int(os.environ.get("TF_KICK_LOSS_BUDGET", "1000"))
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            os.environ.get("TF_KICK_LOSS_OUTPUT", ROOT / "build/kick-loss-comparison")
        ),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.fit and any(args.output.glob("*/search.json")):
        raise ValueError("Existing fits retained: choose a new output directory")
    with ExitStack() as stack:
        renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], "kick-standard", ROOT)
        stack.callback(renderer.close)
        if os.environ.get("TF_KICK_LOSS_RECOVERY") == "1":
            from kick_loss_recovery import run_recovery

            run_recovery(renderer, create_losses, args.output, args.remote_jtfs)
            return
        rate = renderer.sample_rate
        onset = round(renderer.metadata["reference"]["cell"]["onset_seconds"] * rate)
        reference = renderer.reference[onset : onset + round(1.2 * rate)]
        reference = np.pad(reference, (0, round(1.2 * rate) - len(reference)))
        losses = create_losses(reference, rate, args.remote_jtfs, stack)
        starts, bounds = experiment_starts(
            renderer.initial, renderer.metadata["descriptors"]
        )
        manifest = dict(
            metadata=renderer.metadata,
            starts=starts,
            bounds=bounds,
            budget=args.budget,
            duration=1.2,
            publication=False,
            objectives={k: v.specification for k, v in losses.items()},
        )
        save_manifest(args.output, manifest)
        write_wav(args.output / "reference.wav", AudioBuffer(reference, rate))
        audit_losses(losses, reference, rate, args.output)
        score_baselines(renderer, losses, starts, args.output)
        if args.fit:
            run_trials(
                renderer, reference, losses, starts, bounds, args.output, args.budget
            )
        if (args.output / "trials.json").exists():
            from kick_loss_plots import draw_comparison

            draw_comparison(args.output)


if __name__ == "__main__":
    main()
