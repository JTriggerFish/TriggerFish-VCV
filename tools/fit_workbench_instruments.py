"""Fit standard non-kick targets with actual Wasm; no automatic publication.

Run through dev.ps1 fit-instruments. TF_FIT_TARGETS selects comma-separated names,
TF_FIT_ITERATIONS controls iterations per block, TF_FIT_RESUME continues a checked
search, TF_FIT_OUTPUT chooses the experiment root. Audio stays untracked.
"""

import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# Parallel targets each have an independent renderer. Avoid multiplying the
# machine's full BLAS thread pool by the target count.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_rerender import load_checked_start
from triggerfish_percussion.instrument_fit_configuration import configure_fit
from triggerfish_percussion.instrument_fit_stages import metallic_stages, snare_stages
from triggerfish_percussion.metallic_fit_starts import metallic_starts
from triggerfish_percussion.workbench_fit_baseline import original_baseline
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.workbench_multistart import refine_candidate_starts
from triggerfish_percussion.observation_fit_polish import polish_observation

ROOT = Path(__file__).resolve().parents[1]


def fit_target(target, output, iterations, resume):
    output.mkdir(parents=True, exist_ok=True)
    renderer = WorkbenchRenderer(os.environ["EMSDK_NODE"], f"{target}-standard", ROOT)
    try:
        saved = provenance = None
        if (output / "search.json").exists():
            if not resume:
                raise ValueError(
                    "Existing search: set TF_FIT_RESUME=1 or choose a new output"
                )
            saved, provenance = load_checked_start(renderer, output / "search.json")
        seconds, reference, loss = configure_fit(renderer, target, saved, os.environ)
        baseline = original_baseline(renderer, output, seconds)
        write_wav(
            output / "reference.wav", AudioBuffer(reference, renderer.sample_rate)
        )
        baseline_metadata = json.loads((output / "baseline-metadata.json").read_text())
        baseline_parameters = dict(
            zip(
                [item["key"] for item in baseline_metadata["descriptors"]],
                baseline_metadata["values"],
            )
        )
        seed = renderer.metadata["event"]["seed"]
        search = Search(
            renderer,
            loss,
            output,
            seconds,
            target.title(),
            seeds=(
                tuple(saved["training_seeds"])
                if saved
                else (None, (seed + 101) & 0xFFFFFFFF)
            ),
        )
        if saved:
            search.parameters = saved["parameters"]
            search.history = saved["history"]
            search.history.append(
                dict(stage="checked warm start", provenance=provenance)
            )
        print(
            json.dumps(dict(target=target, baseline=loss.diagnostics(baseline))),
            flush=True,
        )
        # Fit one exposed master gain first, then freeze it during shape blocks.
        surface_only = os.environ.get("TF_FIT_PHASE") == "surface"
        if not surface_only:
            search.stage(
                "audible level (visible master only)",
                {"model_level_db": (-36, 0)},
                iterations,
            )
        if target != "snare" and (
            not resume or os.environ.get("TF_FIT_LAYOUTS") == "1"
        ):
            # A raw layout score confounds mode placement with its untuned gain
            # and damping. Compare short refits, not uncalibrated spectra.
            refine_candidate_starts(
                search,
                metallic_starts(search.parameters, reference, renderer.sample_rate),
                lambda _: dict(
                    field_gain=(0.1, 4),
                    body_brightness=(-50, 12),
                    body_decay_seconds_0=(0.1, 30),
                    body_decay_seconds_7=(0.02, 15),
                ),
                count=2,
                iterations=6,
            )
        for name, bounds, regions in (
            snare_stages if target == "snare" else metallic_stages
        )(search.parameters):
            if surface_only and name == "two endpoint damping and energy travel":
                continue
            if name == "modal observation spectrum":
                if os.environ.get("TF_FIT_OBSERVATION") == "autograd":
                    from triggerfish_percussion.torch_observation_fit import (
                        polish_observation_autograd,
                    )

                    polish_observation_autograd(search, iterations)
                else:
                    polish_observation(search, iterations)
            elif bounds:
                search.stage(name, bounds, iterations, regions)
        search.save()
        heldout = []
        for offset in (211, 307, 401):
            fresh = (seed + offset) & 0xFFFFFFFF
            before = loss.diagnostics(
                renderer.render(baseline_parameters, seconds, fresh)
            )
            after = loss.diagnostics(renderer.render(search.parameters, seconds, fresh))
            heldout.append(dict(seed=fresh, baseline=before, candidate=after))
        report = dict(
            target=target,
            baseline=loss.diagnostics(baseline),
            candidate=loss.diagnostics(search.audio(search.parameters)),
            heldout=heldout,
            listening_approved=False,
        )
        (output / "comparison.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(json.dumps(dict(target=target, complete=report)), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    targets = os.environ.get("TF_FIT_TARGETS", "snare,crash,ride,hihat,gong").split(",")
    if not targets or any(
        t not in ("snare", "crash", "ride", "hihat", "gong") for t in targets
    ):
        raise ValueError("Choose non-kick standard targets only")
    if len(set(targets)) != len(targets):
        raise ValueError("Duplicate targets would share an output directory")
    workers = int(os.environ.get("TF_FIT_WORKERS", "1"))
    if not 1 <= workers <= 3:
        raise ValueError("Use one to three independent target workers")

    def run(target):
        return fit_target(
            target,
            ROOT
            / os.environ.get("TF_FIT_OUTPUT", "build/instrument-refits-v1")
            / target,
            int(os.environ.get("TF_FIT_ITERATIONS", "12")),
            os.environ.get("TF_FIT_RESUME") == "1",
        )

    with ThreadPoolExecutor(max_workers=workers) as executor:
        list(executor.map(run, targets))
