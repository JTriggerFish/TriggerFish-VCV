"""Reposition existing handles from reference peaks, without adding DSP."""

import json
import numpy as np
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.modal_fit_initialization import spectral_mode_candidates
from triggerfish_percussion.region_fit_loss import RegionFitLoss
from triggerfish_percussion.workbench_search import Search
from triggerfish_percussion.workbench_fit_baseline import check_reference
from kick_span_refinement import trial_start, assess_trial


def measured_start(renderer, reference, previous):
    proposals = spectral_mode_candidates(reference, renderer.sample_rate)
    proposals = sorted(
        (p for p in proposals if p["frequency"] > 110),
        key=lambda p: p["prominence_db"],
        reverse=True,
    )[:14]
    proposals.sort(key=lambda p: p["frequency"])
    values = dict(previous)
    source = sorted(
        (previous[f"resonance_frequency_{i}"], previous[f"resonance_level_{i}"])
        for i in range(16)
    )
    for index, proposal in enumerate(proposals, 2):
        frequency = proposal["frequency"]
        values.update(
            {
                f"resonance_frequency_{index}": frequency,
                f"resonance_level_{index}": float(
                    np.interp(
                        np.log(frequency),
                        np.log([f for f, l in source]),
                        [l for f, l in source],
                    )
                ),
            }
        )
    return values, proposals


def refine_measured(renderer, audit, reference, output):
    record = json.loads((output / "span-refinement.json").read_text(encoding="utf8"))
    check_reference(record["metadata"], renderer.metadata)
    existing = next(
        (row for row in record["trials"] if row["name"] == "measured-lowpass"), None
    )
    if existing:
        initial, proposals = existing["parameters"], existing["proposals"]
    else:
        previous = next(
            row["parameters"]
            for row in record["trials"]
            if row["name"] == "coverage-lowpass"
        )
        initial, proposals = measured_start(renderer, reference, previous)
    loss = RegionFitLoss(audit, reference, renderer.sample_rate)
    event_seed = renderer.metadata["event"]["seed"]
    seeds = (event_seed, event_seed + 11)
    heldout = tuple(event_seed + i for i in (1, 2, 3))
    directory = output / "measured-lowpass"
    directory.mkdir(parents=True, exist_ok=True)
    search = Search(renderer, loss, directory, 1.2, "measured-lowpass", seeds)
    if existing and (directory / "search.json").exists():
        search.history = json.loads(
            (directory / "search.json").read_text(encoding="utf8")
        )["history"]
    search.parameters = initial
    _, bounds = trial_start(initial, True)
    bounds.update(
        contact_level=(0, 4), thump_level=(1, 4), thump_decay_seconds=(0.15, 0.5)
    )
    search.stage(
        "reference-peak layout; sources and shared damping",
        bounds,
        35,
        difference_step=0.002,
    )
    audio, row = assess_trial(renderer, loss, search.parameters, heldout)
    write_wav(output / "measured-lowpass.wav", AudioBuffer(audio, renderer.sample_rate))
    row.update(
        name="measured-lowpass",
        parameters=search.parameters,
        bounds=bounds,
        proposals=proposals,
        renders=search.evaluations,
        training_seeds=seeds,
        heldout_seeds=heldout,
        search_log=str(directory / "search.json"),
    )
    record["trials"] = [r for r in record["trials"] if r["name"] != row["name"]] + [row]
    (output / "span-refinement.json").write_text(
        json.dumps(record, indent=2), encoding="utf8"
    )
    print(json.dumps(row), flush=True)
