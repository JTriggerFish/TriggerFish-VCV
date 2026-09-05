"""Kick-specific shortlist policy on the reusable exact-renderer multistart."""

from triggerfish_percussion.modal_fit_initialization import extended_modal_starts
from triggerfish_percussion.workbench_multistart import refine_candidate_starts
from kick_fit_stages import stages_for


def refine_modal_alternatives(search, proposals):
    def bounds(values):
        modes = stages_for(values)[1][1]
        return dict(
            modes,
            contact_level=(0, 4),
            thump_level=(0, 4),
            resonance_level=(0, 12),
            resonance_decay_seconds=(0.05, 2),
            resonance_decay_tilt=(-1, 1),
        )

    return refine_candidate_starts(
        search,
        extended_modal_starts(search.parameters, proposals),
        bounds,
        count=3,
        iterations=10,
    )
