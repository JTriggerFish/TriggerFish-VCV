"""Frozen-candidate audit on fresh seeds; never fits or publishes parameters."""

import json
import hashlib
import argparse
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer


def audit_baseline(target, directory, source):
    """Reuse the recorded comparator even after a candidate is published."""
    report_path = directory / "fresh-seed-audit.json"
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding="utf8"))
        if report["target"] != target:
            raise ValueError("Audit target differs from the frozen baseline")
        parameters = report["baseline_parameters"]
    else:
        previous = json.loads(source.read_text(encoding="utf8"))
        parameters = {
            key: value
            for node in previous["instrument"]["nodes"]
            for key, value in node["parameters"].items()
        }
    digest = hashlib.sha256(json.dumps(parameters, sort_keys=True).encode()).hexdigest()
    if report_path.exists() and report.get("baseline_sha256", digest) != digest:
        raise ValueError("Frozen audit baseline hash mismatch")
    return parameters, digest


def verify(target, directory, offsets=(1009, 2017, 3011)):
    if len(set(offsets)) != len(offsets) or any(not 0 < x < 2**32 for x in offsets):
        raise ValueError("Audit offsets must be distinct nonzero 32-bit integers")
    root = Path.cwd()
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], target + "-standard", root)
    try:
        saved = verify_candidate(renderer, directory)
        # The previous workbench controls are the comparator. The corrected
        # source onset applies to BOTH waveforms; don't compare obsolete scores.
        baseline, baseline_hash = audit_baseline(
            target, directory, root / "workbench/web" / f"{target}_calibration.fit.json"
        )
        seconds = saved["duration_seconds"]
        reference = aligned_reference(renderer, seconds)
        loss = MetallicBalanceLoss(
            reference, renderer.sample_rate, contrast_weighting="erb", fast_attack=True
        )
        seed = renderer.metadata["event"]["seed"]
        rows = []
        for offset in (0, *offsets):
            actual = (seed + offset) & 0xFFFFFFFF
            before = renderer.render(baseline, seconds, actual)
            after = renderer.render(saved["parameters"], seconds, actual)
            rows.append(
                dict(
                    seed=actual,
                    role="standard" if offset == 0 else "audit",
                    baseline=loss.diagnostics(before),
                    candidate=loss.diagnostics(after),
                    baseline_peak=float(np.max(np.abs(before))),
                    candidate_peak=float(np.max(np.abs(after))),
                )
            )
        report = dict(
            target=target,
            objective=loss.specification,
            units=loss.units,
            renderer=renderer.metadata,
            baseline_parameters=baseline,
            baseline_sha256=baseline_hash,
            seed_offsets=offsets,
            candidate_parameters=saved["parameters"],
            seeds=rows,
            listening_approved=False,
        )
        (directory / "fresh-seed-audit.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(json.dumps(dict(target=target, seeds=rows)), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("crash", "ride", "gong", "hihat"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--seed-offsets", type=int, nargs=3, default=(1009, 2017, 3011))
    arguments = parser.parse_args()
    verify(arguments.target, arguments.directory, tuple(arguments.seed_offsets))
