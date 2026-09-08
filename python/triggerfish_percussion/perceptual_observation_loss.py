"""Optional auraloss autograd adapter for exact-render observation fitting."""

import copy
import numpy as np
import torch


class TorchAuralossMel:
    """Differentiate the library objective, not a replacement mel implementation.

    Double precision improves finite-difference validation of amplitude gradients.
    Keep a separate module so the public float32 scoring path remains unchanged.
    """

    def __init__(self, source):
        self.source = source
        self.loss = copy.deepcopy(source.loss).double()
        self.target = source.target.double()

    def __call__(self, samples):
        return self.loss(samples[None, None], self.target)

    def validate(self, audio):
        expected = self.source.score(audio)
        actual = float(self(torch.tensor(audio, dtype=torch.float64)))
        relative = abs(actual - expected) / max(expected, 1)
        if not np.isfinite((expected, actual, relative)).all() or relative > 1e-5:
            raise ValueError(f"Differentiable auraloss differs: {actual=}, {expected=}")
        return dict(expected=expected, actual=actual, relative_error=relative)
