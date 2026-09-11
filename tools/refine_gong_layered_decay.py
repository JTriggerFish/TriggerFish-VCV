"""Small exact-render shared-damping screen after grouped observation fitting."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch

from triggerfish_percussion.fit_reference import aligned_reference
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from fit_stretched_gong import checkpoint
from refine_gong_layered_bloom import LayeredBloomLoss


def run(args):
    torch.set_num_threads(1)
    r = WorkbenchRenderer(os.environ["EMSDK_NODE"], "gong-standard", Path.cwd())
    try:
        source = json.loads((args.source / "search.json").read_text())
        base = source["parameters"]
        user = json.loads(args.user.read_text())["parameters"]
        reference = aligned_reference(r, 6)
        loss = LayeredBloomLoss(reference, r.render(user, 6), r.sample_rate)
        trials = [dict(base)]
        for high_t60 in (2.14, 3.5, 5, 7):
            for rate in (3, 4, 6):
                for concentration in (0.2, 0.35):
                    trials.append(
                        dict(
                            base,
                            body_decay_seconds_7=high_t60,
                            bloom_rate=rate,
                            bloom_energy_acceleration=concentration,
                        )
                    )
        args.output.mkdir(parents=True, exist_ok=True)
        rows = []
        for i, parameters in enumerate(trials):
            audio = [r.render(parameters, 6, s) for s in (1675, 1982)]
            rows.append(
                dict(
                    parameters=parameters,
                    score=float(np.mean([loss.score(a) for a in audio])),
                    diagnostics=[loss.diagnostics(a) for a in audio],
                )
            )
            (args.output / "screen.json").write_text(json.dumps(rows, indent=2))
            if i % 6 == 0:
                print(
                    json.dumps(dict(trial=i, best=min(x["score"] for x in rows))),
                    flush=True,
                )
        best = min(rows, key=lambda x: x["score"])
        search = checkpoint(
            r,
            loss,
            args.output,
            "Gong — tuned body and bloom",
            best["parameters"],
            reference,
            source["history"]
            + [
                dict(
                    stage="shared high damping screen",
                    trials=len(rows),
                    before=rows[0]["score"],
                    after=best["score"],
                    objective=loss.specification,
                    fitting_seeds=[1675, 1982],
                )
            ],
        )
        verify_candidate(r, args.output)
        audit = [
            dict(
                seed=s, diagnostics=loss.diagnostics(r.render(best["parameters"], 6, s))
            )
            for s in (1675, 1982, 2586, 3276)
        ]
        (args.output / "audit.json").write_text(json.dumps(audit, indent=2))
        print(
            json.dumps(dict(before=rows[0]["score"], after=best["score"], audit=audit)),
            flush=True,
        )
    finally:
        r.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("user", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
