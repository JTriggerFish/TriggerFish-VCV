"""Listening-directed short-contact alternatives, using the unchanged C++ voice.

Noise <= 70 ms is an explicit experimental prior from user feedback, not a DSP
limit or a measured source separation. Pitch remains free, including sub-bass.
"""

import json

from triggerfish_percussion.workbench_search import Search
from kick_fit_stages import stages_for


def refine_attack(search, impact_starts=False):
    initial = dict(search.parameters)
    bounds = {
        key: value for _, stage in stages_for(initial) for key, value in stage.items()
    }
    bounds["contact_noise_decay_seconds"] = (0.005, 0.07)
    # Fit each constrained start independently: it need not beat the old,
    # unconstrained patch before it has a chance to rebalance the sources.
    candidates = []
    starts = [
        dict(thump_pitch_hz=pitch, contact_noise_decay_seconds=decay)
        for pitch, decay in ((27, 0.035), (40, 0.035), (55, 0.02))
    ]
    prefix = "short-contact"
    if impact_starts:
        prefix = "short-impact"
        starts = [
            dict(
                contact_width_seconds=width,
                contact_noise_level=0.3,
                contact_noise_decay_seconds=0.025,
            )
            for width in (0.001, 0.003)
        ]
    for index, overrides in enumerate(starts):
        directory = search.output / f"{prefix}-{index}"
        directory.mkdir(parents=True, exist_ok=True)
        child = Search(
            search.renderer,
            search.loss,
            directory,
            search.seconds,
            search.name,
            search.seeds,
        )
        child.parameters = dict(initial, **overrides)
        source_bounds = {
            key: value
            for key, value in bounds.items()
            if not key.startswith("resonance_frequency_")
            and not key.startswith("resonance_level_")
        }
        child.stage(
            "short contact: rebalance sources and shared damping", source_bounds, 22
        )
        child.stage("short contact: joint sources and existing modes", bounds, 30)
        candidates.append((f"{prefix}-{index}", dict(child.parameters)))
    search.screen_candidates(
        "short-contact refits against unchanged baseline", candidates
    )
    if search.parameters != initial:
        search.stage(
            "short contact: fine refinement", bounds, 30, difference_step=0.001
        )
    report = dict(
        prior="Contact noise base T60 <= 70 ms; listening-directed experiment",
        before=search.loss.diagnostics(search.audio(initial)),
        after=search.loss.diagnostics(search.audio(search.parameters)),
        candidates=[dict(name=name, parameters=values) for name, values in candidates],
    )
    (search.output / "attack-refinement.json").write_text(json.dumps(report, indent=2))
