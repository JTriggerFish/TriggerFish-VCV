"""Exact-renderer recovery test and saved-fit validation, without real-data claims."""

import json

import numpy as np

from triggerfish_percussion.short_drum_fit_loss import ShortDrumLoss
from triggerfish_percussion.workbench_search import Search


def self_test(renderer, output):
    """Recover known C++ membrane pitch and damping before blaming the reference."""
    output.mkdir(parents=True, exist_ok=True)
    known = dict(
        renderer.initial,
        resonance_frequency_0=49,
        resonance_decay_seconds=0.7,
        resonance_level=2,
        contact_level=0.05,
        thump_level=0,
        tension_octaves=0,
    )
    known.update({f"resonance_level_{i}": -72 for i in range(1, 16)})
    reference = renderer.render(known, 1.2)
    loss = ShortDrumLoss(reference, renderer.sample_rate)
    search = Search(
        renderer, loss, output, seconds=1.2, name="Synthetic recovery fixture"
    )
    search.parameters = dict(
        known, resonance_frequency_0=43, resonance_decay_seconds=0.95
    )
    search.stage(
        "known C++ pitch and damping recovery",
        {"resonance_frequency_0": (35, 65), "resonance_decay_seconds": (0.3, 1.2)},
        30,
    )
    recovered = search.parameters
    if (
        abs(recovered["resonance_frequency_0"] - 49) > 0.15
        or abs(recovered["resonance_decay_seconds"] - 0.7) > 0.01
    ):
        raise RuntimeError("Exact C++ synthetic parameter recovery failed")
    print(
        json.dumps(
            dict(
                self_test="passed",
                pitch=recovered["resonance_frequency_0"],
                t60=recovered["resonance_decay_seconds"],
            )
        ),
        flush=True,
    )
    recover_contact_noise(renderer, output / "contact-noise")


def recover_contact_noise(renderer, output):
    """Recover the actual recipe's coloured-noise level and finite fade time."""
    output.mkdir(parents=True, exist_ok=True)
    known = dict(
        renderer.initial,
        resonance_level=0,
        thump_level=0,
        contact_level=0.5,
        contact_noise_level=1.2,
        contact_noise_decay_seconds=0.18,
        contact_colour=0.7,
    )
    reference = renderer.render(known, 1.2)
    search = Search(
        renderer,
        ShortDrumLoss(reference, renderer.sample_rate),
        output,
        seconds=1.2,
        name="Contact noise recovery fixture",
    )
    search.parameters = dict(
        known, contact_noise_level=0.5, contact_noise_decay_seconds=0.07
    )
    search.stage(
        "known C++ contact noise recovery",
        {"contact_noise_level": (0.1, 3), "contact_noise_decay_seconds": (0.03, 0.4)},
        30,
    )
    recovered = search.parameters
    if (
        abs(recovered["contact_noise_level"] - 1.2) > 0.03
        or abs(recovered["contact_noise_decay_seconds"] - 0.18) > 0.003
    ):
        raise RuntimeError("Exact C++ contact-noise recovery failed")
    print(json.dumps(dict(noise_self_test="passed", parameters=recovered)), flush=True)


def validate_candidate(renderer, search, loss, output, seconds):
    fit = renderer.request(command="snapshot", parameters=search.parameters)["fit"]
    restored = renderer.decode(
        renderer.request(command="renderSnapshot", fit=fit, seconds=seconds)["pcm"]
    )
    if not np.array_equal(restored, search.audio(search.parameters)):
        raise RuntimeError("Fit JSON does not reproduce the search render")
    heldout = [
        loss.diagnostics(
            renderer.render(
                search.parameters, seconds, renderer.metadata["event"]["seed"] + i
            )
        )
        for i in (1, 2, 3)
    ]
    (output / "heldout-seeds.json").write_text(json.dumps(heldout, indent=2))
    return restored
