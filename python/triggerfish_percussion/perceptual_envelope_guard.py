"""Reference-fixed band/decay and attack constraints for observation fitting."""

import numpy as np
import torch
from .band_decay_shape_loss import BandDecayShapeLoss
from .metallic_balance_features import ATTACK_BINS
from .torch_metallic_features import TorchMetallicFeatures


class PerceptualEnvelopeGuard:
    """Protect each seed separately, without normalizing either output waveform.

    Band absolute and relative-decay errors cover the reference's usable decay
    mask (from 0.2 s). Five additional raw-power bins protect the initial 100 ms.
    These are non-regression constraints, not a claim of perceptual equivalence.
    """

    def __init__(self, basis, reference, rate, baseline_audio, tolerance_db=0.05):
        if not np.isfinite(tolerance_db) or tolerance_db < 0:
            raise ValueError("Guard tolerance must be finite and nonnegative")
        self.basis = basis
        self.shape = BandDecayShapeLoss(reference, rate)
        self.reference_db = torch.tensor(
            10 * np.log10(np.maximum(self.shape.power(reference), self.shape.floor))
        )
        self.reference_relative = torch.tensor(self.shape.target)
        self.anchor = torch.tensor(self.shape.anchor)
        self.sections = [
            slice(round(a * rate), round(b * rate)) for a, b in ATTACK_BINS
        ]
        self.attack_floor = max(float(np.max(reference**2)) * 1e-10, 1e-20)
        self.attack_db = torch.tensor(
            [
                10 * np.log10(max(float(np.mean(reference[s] ** 2)), self.attack_floor))
                for s in self.sections
            ]
        )
        self.prepared = [
            (torch.tensor(audio), torch.tensor(columns))
            for audio, columns in basis.bases.values()
        ]
        self.initial = torch.tensor(basis.amplitudes)
        if set(baseline_audio) != set(basis.bases):
            raise ValueError("Guard comparator must cover the exact training seeds")
        self.limits = torch.cat(
            [
                (
                    self.errors(torch.tensor(baseline_audio[seed])).sqrt()
                    + tolerance_db
                ).square()
                for seed in basis.bases
            ]
        ).detach()
        self.cached = None
        self.specification = dict(
            version="perceptual-envelope-guard-v1",
            tolerance_db=tolerance_db,
            decay=self.shape.specification,
            attack_bins_seconds=ATTACK_BINS,
            protects="per-band absolute and decay RMS; per-bin attack power",
        )

    def errors(self, samples):
        if (
            samples.shape != (self.shape.frames,)
            or samples.is_complex()
            or not torch.isfinite(samples).all()
        ):
            raise ValueError("Guard audio must be finite mono and match the reference")
        samples = samples.to(dtype=torch.float64)
        transformed = TorchMetallicFeatures.stft(samples, 4096, 512).abs().square()
        power = torch.stack(
            [transformed[mask].sum(dim=0) for mask in self.shape.bin_masks]
        )
        db = 10 * torch.log10(power.clamp_min(self.shape.floor))
        relative = db - db[:, self.anchor].mean(dim=1, keepdim=True)
        rows = []
        for index, mask in enumerate(self.shape.mask):
            if mask.sum() < 4:
                continue
            rows.extend(
                (
                    (db[index, mask] - self.reference_db[index, mask]).square().mean(),
                    (relative[index, mask] - self.reference_relative[index, mask])
                    .square()
                    .mean(),
                )
            )
        attack = torch.stack(
            [
                10
                * torch.log10(samples[s].square().mean().clamp_min(self.attack_floor))
                for s in self.sections
            ]
        )
        return torch.cat((torch.stack(rows), (attack - self.attack_db).square()))

    def margins(self, amplitudes):
        return self.limits - torch.cat(
            [
                self.errors(audio + (amplitudes - self.initial) @ columns)
                for audio, columns in self.prepared
            ]
        )

    def values(self, amplitudes):
        with torch.no_grad():
            return self.margins(torch.tensor(amplitudes, dtype=torch.float64)).numpy()

    def evaluate(self, amplitudes):
        token = np.asarray(amplitudes).tobytes()
        if self.cached is not None and self.cached[0] == token:
            return self.cached[1:]
        weights = torch.tensor(amplitudes, dtype=torch.float64, requires_grad=True)
        margins = self.margins(weights)
        jacobian = torch.stack(
            [
                torch.autograd.grad(value, weights, retain_graph=True)[0]
                for value in margins
            ]
        )
        result = (margins.detach().numpy(), jacobian.numpy())
        if not all(np.isfinite(value).all() for value in result):
            raise ValueError("Nonfinite envelope guard")
        self.cached = (token, *result)
        return result
