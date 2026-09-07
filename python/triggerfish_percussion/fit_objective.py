"""Preserve a recorded objective when resuming or polishing a search."""

import json

from .metallic_balance_loss import MetallicBalanceLoss
from .metallic_fit_loss import MetallicFitLoss


def check_objective(saved, loss):
    """Changing measurements requires a new experiment, not overwritten history."""
    expected = json.dumps(saved.get("objective_specification"), sort_keys=True)
    actual = json.dumps(getattr(loss, "specification", None), sort_keys=True)
    if saved.get("objective") != type(loss).__name__ or expected != actual:
        raise ValueError("Fitting objective changed; use a new experiment directory")


def saved_metallic_loss(saved, reference, rate):
    """Recreate only known, fully specified metallic measurements."""
    name = saved.get("objective")
    if name == "MetallicBalanceLoss":
        specification = saved.get("objective_specification") or {}
        loss = MetallicBalanceLoss(
            reference,
            rate,
            contrast_weighting=specification.get("contrast_weighting", "linear"),
            fast_attack=specification.get("fast_attack", False),
        )
    elif name == "MetallicFitLoss":
        loss = MetallicFitLoss(reference, rate)
    else:
        raise ValueError(f"Unsupported saved metallic objective: {name}")
    check_objective(saved, loss)
    return loss
