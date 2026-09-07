"""Test coordinated control reachability, never automatically publish a preset."""

import json
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.region_fit_loss import RegionFitLoss
from triggerfish_percussion.workbench_search import Search
from kick_control_span import coverage_layout


def trial_start(baseline, filtered):
    values = coverage_layout(baseline)
    values.update(
        contact_width_seconds=0.003,
        contact_noise_decay_seconds=0.05,
        contact_level=0.1,
        resonance_decay_seconds=0.6,
    )
    bounds = dict(
        contact_level=(0, 1),
        contact_width_seconds=(0.0005, 0.025),
        contact_noise_level=(0, 4),
        contact_noise_decay_seconds=(0.005, 0.4),
        contact_colour=(0, 1),
        resonance_level=(0, 12),
        resonance_decay_seconds=(0.05, 1),
        resonance_decay_tilt=(-0.5, 1),
    )
    bounds.update({f"resonance_level_{i}": (-60, 6) for i in range(16) if i != 1})
    if filtered:
        values.update(
            equalizer_mode=1, low_cut_hz=10, high_cut_hz=2500, colour_gain_db=0
        )
        bounds["high_cut_hz"] = (600, 8000)
    return values, bounds


def assess_trial(renderer, loss, parameters, heldout_seeds):
    audio = renderer.render(parameters, 1.2)
    row = loss.diagnostics(audio)
    row["heldout"] = [
        loss.diagnostics(renderer.render(parameters, 1.2, seed))
        for seed in heldout_seeds
    ]
    return audio, row


def refine_span(renderer, audit, reference, output):
    loss = RegionFitLoss(audit, reference, renderer.sample_rate)
    event_seed = renderer.metadata["event"]["seed"]
    training_seeds = (event_seed, event_seed + 11)
    heldout_seeds = tuple(event_seed + i for i in (1, 2, 3))
    baseline = loss.diagnostics(renderer.render(renderer.initial, 1.2))
    results = []
    for filtered in (False, True):
        label = "coverage-lowpass" if filtered else "coverage-bypass"
        directory = output / label
        directory.mkdir(parents=True, exist_ok=True)
        search = Search(renderer, loss, directory, 1.2, label, training_seeds)
        search.parameters, bounds = trial_start(renderer.initial, filtered)
        search.stage(
            "coordinated control-span diagnostic", bounds, 35, difference_step=0.002
        )
        audio, row = assess_trial(renderer, loss, search.parameters, heldout_seeds)
        write_wav(output / f"{label}.wav", AudioBuffer(audio, renderer.sample_rate))
        row.update(
            name=label,
            parameters=search.parameters,
            bounds=bounds,
            renders=search.evaluations,
            training_seeds=training_seeds,
            heldout_seeds=heldout_seeds,
            search_log=str(directory / "search.json"),
        )
        results.append(row)
        report = dict(
            metadata=renderer.metadata,
            baseline_band=baseline["band"],
            baseline_shape=baseline["shape"],
            trials=results,
        )
        (output / "span-refinement.json").write_text(
            json.dumps(report, indent=2), encoding="utf8"
        )
        print(json.dumps(row), flush=True)
