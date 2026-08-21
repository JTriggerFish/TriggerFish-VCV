"""Differentiable reference models for TriggerFish reverb research.

The recovered VFM model is an offline research tool. It is not imported or
executed by the Rack plugin.
"""

from .objectives import aggregate_resonance, resonance_loss, resonance_per_control
from .velvet import DifferentiableVelvetReverb, VelvetControls

__all__ = [
    "DifferentiableVelvetReverb",
    "VelvetControls",
    "aggregate_resonance",
    "resonance_loss",
    "resonance_per_control",
]
