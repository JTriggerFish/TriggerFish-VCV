"""Explicitly migrate a parameter checkpoint to the current renderer.

Old WAVs/scores are NOT carried forward after a DSP change. Validate the old
snapshot's internal consistency, record its renderer hash, then render, score
and verify a new complete checkpoint. No preset publication is automatic.
"""

import argparse
import json
import os
from pathlib import Path

import torch
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import check_saved_snapshot, verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.torch_observation_fit import polish_observation_autograd
from triggerfish_percussion.workbench_fit_baseline import check_reference
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search


def run(args):
    torch.set_num_threads(1)
    source = json.loads((args.source / "search.json").read_text(encoding="utf8"))
    snapshot = json.loads(
        (args.source / "candidate.fit.json").read_text(encoding="utf8")
    )
    check_saved_snapshot(snapshot, source)
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], f"{args.target}-standard", Path.cwd()
    )
    try:
        check_reference(source["metadata"], renderer.metadata)
        overrides = dict(item.split("=", 1) for item in args.set)
        if (set(source["parameters"]) | set(overrides)) != set(renderer.initial):
            raise ValueError(
                "Explicit --set required for each new control; unknown controls rejected"
            )
        seconds = source["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        loss = MetallicBalanceLoss(reference, renderer.sample_rate, "erb", True)
        args.output.mkdir(parents=True, exist_ok=True)
        search = Search(renderer, loss, args.output, seconds, args.name, (None,))
        search.parameters = dict(source["parameters"])
        for key, value in overrides.items():
            search.parameters[key] = float(value)
        search.history.append(
            dict(
                parent=str(args.source.resolve()),
                overrides=overrides,
                previous_renderer=source["metadata"]["rendererSha256"],
                current_renderer=renderer.metadata["rendererSha256"],
                old_scores_reused=False,
            )
        )
        search.audio(search.parameters)  # Validates every value before optimization.
        if args.polish:
            polish_observation_autograd(search, iterations=args.polish)
        search.save()
        write_wav(
            args.output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        verify_candidate(renderer, args.output)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=["crash", "gong"])
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--set", action="append", default=[])
    parser.add_argument("--polish", type=int, default=0)
    run(parser.parse_args())
