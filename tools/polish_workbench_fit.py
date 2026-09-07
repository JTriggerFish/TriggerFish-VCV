"""Post-topology observation refinement; actual C++ checks before and after."""

import json
import os
from pathlib import Path
from itertools import product

from triggerfish_percussion.audio_io import read_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_objective import saved_metallic_loss
from triggerfish_percussion.workbench_fit_baseline import check_reference
from triggerfish_percussion.observation_fit_polish import polish_observation
from triggerfish_percussion.metallic_fit_alternative import (
    late_ridge_parameters,
    refit_alternative,
)
from triggerfish_percussion.workbench_renderer import WorkbenchRenderer
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.fit_rerender import load_checked_start


def polish(target, directory):
    renderer = WorkbenchRenderer(
        os.environ["EMSDK_NODE"], target + "-standard", Path.cwd()
    )
    try:
        saved = verify_candidate(renderer, directory)
        reference = read_wav(directory / "reference.wav").mono().samples
        baseline = read_wav(directory / "baseline.wav").mono().samples
        metadata = json.loads((directory / "baseline-metadata.json").read_text())
        check_reference(metadata, renderer.metadata)
        if len(baseline) != round(saved["duration_seconds"] * renderer.sample_rate):
            raise ValueError("Baseline duration differs from the saved fit")
        baseline_parameters = dict(
            zip([d["key"] for d in metadata["descriptors"]], metadata["values"])
        )
        loss = saved_metallic_loss(saved, reference, renderer.sample_rate)
        search = Search(
            renderer,
            loss,
            directory,
            saved["duration_seconds"],
            target.title(),
            tuple(saved["training_seeds"]),
        )
        search.parameters, search.history = saved["parameters"], saved["history"]
        if os.environ.get("TF_FIT_RIDGE_DAMPING") == "1":
            coherent, _ = load_checked_start(
                renderer, directory / "coherent-tail-layout" / "search.json"
            )
            values = dict(
                coherent["parameters"],
                body_decay_active_1=1,
                body_decay_frequency_1=1000.0,
                body_decay_seconds_1=24.0,
                body_decay_seconds_0=30.0,
                body_decay_seconds_7=0.1,
            )
            refit_alternative(
                search,
                values,
                "coherent-tail-damping",
                int(os.environ.get("TF_FIT_ITERATIONS", "8")),
            )
        if os.environ.get("TF_FIT_LATE_RIDGES") == "1":
            refit_alternative(
                search,
                late_ridge_parameters(
                    search.parameters,
                    reference,
                    renderer.sample_rate,
                    low=60 if target == "crash" else 250,
                ),
                "coherent-tail-layout",
                int(os.environ.get("TF_FIT_ITERATIONS", "8")),
            )
        if os.environ.get("TF_FIT_LOW_FED_ALTERNATIVE") == "1":
            refit_alternative(
                search,
                dict(
                    search.parameters,
                    body_excitation_centre=200.0,
                    body_brightness=-24.0,
                    bloom_rate=8.0,
                    body_decay_seconds_7=2.0,
                ),
                "low-fed-bloom-layout",
                int(os.environ.get("TF_FIT_ITERATIONS", "8")),
            )
        if os.environ.get("TF_FIT_ENERGY_SCREEN") == "1":
            starts = []
            for centre, tilt, rate, high_decay in product(
                (200.0, 500.0),
                (-12.0, -24.0),
                (2.0, 4.0, 8.0),
                (
                    search.parameters["body_decay_seconds_7"],
                    search.parameters["body_decay_seconds_7"] * 0.5,
                ),
            ):
                values = dict(
                    search.parameters,
                    body_excitation_centre=centre,
                    body_brightness=tilt,
                    bloom_rate=rate,
                    body_decay_seconds_7=max(0.02, high_decay),
                )
                starts.append(
                    (
                        f"energy centre {centre}, tilt {tilt}, rate {rate}, high T60 {high_decay}",
                        values,
                    )
                )
            search.screen_candidates("low-fed energy travel regimes", starts)
        if os.environ.get("TF_FIT_TEXTURE_SCREEN") == "1":
            # Jointly cross the coherent/diffuse boundary. Independent tiny
            # perturbations cannot find a clean bank if other diffusers remain on.
            starts = []
            for (phase, exchange, transfer), spread in product(
                ((0, 0, 0), (0.01, 0.02, 0.1), (0.05, 0.05, 0.3)), (1.0, 3.0, 6.0)
            ):
                values = dict(
                    search.parameters,
                    field_phase_bandwidth=phase,
                    field_exchange=exchange,
                    bloom_phase_diffusion=transfer,
                    field_packet_spread=spread,
                )
                starts.append(
                    (
                        f"phase {phase}, exchange {exchange}, transfer {transfer}, spread {spread}",
                        values,
                    )
                )
            search.screen_candidates("joint coherence regimes", starts)
        print(
            json.dumps(
                dict(target=target, phase="preparing validated observation basis")
            ),
            flush=True,
        )
        if not (
            os.environ.get("TF_FIT_RIDGE_DAMPING") == "1"
            or os.environ.get("TF_FIT_LATE_RIDGES") == "1"
            or os.environ.get("TF_FIT_LOW_FED_ALTERNATIVE") == "1"
        ):
            polish_observation(search, int(os.environ.get("TF_FIT_ITERATIONS", "20")))
        heldout = []
        for offset in (211, 307, 401):
            seed = (renderer.metadata["event"]["seed"] + offset) & 0xFFFFFFFF
            heldout.append(
                dict(
                    seed=seed,
                    baseline=loss.diagnostics(
                        renderer.render(baseline_parameters, search.seconds, seed)
                    ),
                    candidate=loss.diagnostics(
                        renderer.render(search.parameters, search.seconds, seed)
                    ),
                )
            )
        report = dict(
            target=target,
            baseline=loss.diagnostics(baseline),
            candidate=loss.diagnostics(search.audio(search.parameters)),
            heldout=heldout,
            listening_approved=False,
        )
        (directory / "comparison.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        verify_candidate(renderer, directory)
        print(json.dumps(dict(target=target, complete=report)), flush=True)
    finally:
        renderer.close()


if __name__ == "__main__":
    targets = os.environ.get("TF_FIT_TARGETS", "crash,ride,hihat,gong").split(",")
    if any(target not in ("crash", "ride", "hihat", "gong") for target in targets):
        raise ValueError("Observation polishing supports metallic targets only")
    root = Path(os.environ.get("TF_FIT_OUTPUT", "build/instrument-refits-v1"))
    for target in targets:
        polish(target, root / target)
