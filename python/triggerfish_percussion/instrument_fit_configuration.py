"""Explicit experiment settings, with saved settings winning over defaults."""

import numpy as np

from .drum_balance_loss import DrumBalanceLoss
from .fit_objective import check_objective
from .fit_reference import aligned_reference
from .metallic_balance_loss import MetallicBalanceLoss
from .metallic_fit_loss import MetallicFitLoss


def configure_fit(renderer, target, saved, environment):
    """Return duration, aligned reference and loss without writing any files."""
    previous = (saved or {}).get("objective_specification") or {}
    default = (
        "trajectory-v1"
        if (saved or {}).get("objective") == "MetallicFitLoss"
        else "balance-v2"
    )
    objective = environment.get("TF_FIT_OBJECTIVE", default)
    if objective not in ("trajectory-v1", "balance-v2"):
        raise ValueError("Choose trajectory-v1 or balance-v2 objective")
    seconds = (
        1.2
        if target == "snare"
        else float(
            environment.get("TF_FIT_SECONDS", (saved or {}).get("duration_seconds", 6))
        )
    )
    if not np.isfinite(seconds) or (target != "snare" and not 3.1 <= seconds <= 30):
        raise ValueError("Metallic fitting duration must be between 3.1 and 30 seconds")
    if saved and seconds != saved["duration_seconds"]:
        raise ValueError("Fitting duration changed; use a new experiment directory")
    reference = aligned_reference(renderer, seconds)
    if target == "snare":
        loss = DrumBalanceLoss(reference, renderer.sample_rate)
    elif objective == "trajectory-v1":
        loss = MetallicFitLoss(reference, renderer.sample_rate)
    else:
        fast = environment.get(
            "TF_FIT_FAST_ATTACK", "1" if previous.get("fast_attack", not saved) else "0"
        )
        if fast not in ("0", "1"):
            raise ValueError("TF_FIT_FAST_ATTACK must be 0 or 1")
        loss = MetallicBalanceLoss(
            reference,
            renderer.sample_rate,
            contrast_weighting=environment.get(
                "TF_FIT_CONTRAST", previous.get("contrast_weighting", "erb")
            ),
            fast_attack=fast == "1",
        )
    if saved:
        check_objective(saved, loss)
    return seconds, reference, loss
