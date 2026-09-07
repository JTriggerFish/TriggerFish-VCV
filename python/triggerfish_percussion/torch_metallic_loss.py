"""Autograd equivalent of MetallicBalanceLoss, with identical fixed weights."""

import numpy as np
import torch

from .torch_metallic_features import TorchMetallicFeatures


class TorchMetallicLoss:
    def __init__(self, loss):
        self.source = loss
        self.features = TorchMetallicFeatures(loss.features)
        self.target = {key: torch.tensor(value) for key, value in loss.target.items()}
        self.db = {key: torch.tensor(value) for key, value in loss.db.items()}
        self.weights = {key: torch.tensor(value) for key, value in loss.weights.items()}
        self.contrast = torch.tensor(loss.contrast)

    @staticmethod
    def mean_square(error, weight):
        return (error.square() * weight).sum() / weight.sum()

    @staticmethod
    def safe_sqrt(value):
        """Exact zero magnitude, with a finite zero derivative at silence."""
        root = value.clamp_min(torch.finfo(value.dtype).tiny).sqrt()
        return torch.where(value > 0, root, torch.zeros_like(root))

    def __call__(self, samples):
        values = self.features(samples)
        db = {
            key: 10 * torch.log10(value.clamp_min(self.source.floors[key]))
            for key, value in values.items()
        }
        contrast = db["spectrum"] - self.features.smooth_frequency(db["spectrum"], 8)
        parts = dict(envelope=0.0, linear_spectrum=0.0, contrast=0.0)
        for r, (a, b) in enumerate(self.source.features.regions):
            times = self.source.features.times
            masks = [(times >= a) & (times < b)]
            if r == 0:
                masks = [masks[0] & (times < 0.03), masks[0] & (times >= 0.03)]
            for mask in masks:
                parts["envelope"] += self.mean_square(
                    db["envelope"][:, mask] - self.db["envelope"][:, mask],
                    self.weights["envelope"][:, mask],
                ) / (5 * len(masks))
            ref = self.target["spectrum"][r].sqrt()
            denominator = max(
                float(ref.square().sum()), self.source.floors["spectrum"] * len(ref)
            )
            parts["linear_spectrum"] += (
                400
                * (self.safe_sqrt(values["spectrum"][r]) - ref).square().sum()
                / denominator
                / 5
            )
            parts["contrast"] += (
                self.mean_square(
                    contrast[r] - self.contrast[r], self.weights["spectrum"][r]
                )
                / 5
            )
        parts["attack"] = self.mean_square(
            db["attack"] - self.db["attack"], self.weights["attack"]
        )
        if self.source.specification["fast_attack"]:
            parts["attack"] = 0.5 * (
                parts["attack"]
                + (db["transient"] - self.db["transient"]).square().mean()
            )
        return sum(
            self.source.specification["shares"][key] * value
            for key, value in parts.items()
        )

    def validate(self, audio):
        expected = float(np.linalg.norm(self.source.residual(audio)) ** 2)
        actual = float(self(torch.tensor(audio)))
        relative = abs(actual - expected) / max(expected, 1)
        if not np.isfinite((expected, actual, relative)).all() or relative > 1e-5:
            raise ValueError(
                f"Differentiable measurement differs: {actual=}, {expected=}"
            )
        return dict(expected=expected, actual=actual, relative_error=relative)
