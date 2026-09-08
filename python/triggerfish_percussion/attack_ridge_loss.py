"""Optional library MR-STFT loss with linear frequency bins for attack ridges."""

from importlib.metadata import version

import numpy as np

from .perceptual_fit_losses import ScalarAudioLoss


class AttackRidgeLoss(ScalarAudioLoss):
    def __init__(self, reference, rate, seconds=0.4):
        import torch
        from auraloss.freq import MultiResolutionSTFTLoss

        reference = np.asarray(reference)
        if not np.isfinite((rate, seconds)).all() or rate <= 0 or seconds <= 0:
            raise ValueError("Invalid attack sample rate or duration")
        self.frames = round(rate * seconds)
        if self.frames <= 8192 or len(reference) < self.frames:
            raise ValueError("Attack window must contain more than 8192 samples")
        if reference.ndim != 1 or not np.isfinite(reference[: self.frames]).all():
            raise ValueError("Invalid mono attack reference")
        self.torch = torch
        self.target = torch.tensor(reference[: self.frames], dtype=torch.float32)[
            None, None
        ]
        self.settings = dict(
            fft_sizes=[4096, 8192, 16384],
            hop_sizes=[1024, 2048, 4096],
            win_lengths=[4096, 8192, 16384],
            scale=None,
            scale_invariance=False,
            w_sc=1.0,
            w_log_mag=1.0,
            eps=1e-10,
        )
        self.loss = MultiResolutionSTFTLoss(**self.settings)
        self.specification = dict(
            implementation="auraloss",
            version=version("auraloss"),
            attack_seconds=seconds,
            sample_rate=rate,
            settings=self.settings,
        )

    def score(self, samples):
        samples = np.asarray(samples)
        if (
            samples.ndim != 1
            or len(samples) < self.frames
            or not np.isfinite(samples[: self.frames]).all()
        ):
            raise ValueError("Invalid attack samples")
        with self.torch.inference_mode():
            candidate = self.torch.tensor(
                samples[: self.frames], dtype=self.torch.float32
            )[None, None]
            return float(self.loss(candidate, self.target))
