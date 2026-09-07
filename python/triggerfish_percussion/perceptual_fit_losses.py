"""Optional published audio losses; no dependency from the Rack/runtime DSP.

Each callable returns its native scalar objective, without level matching.
The residual adapter is for diagnostics only, not least-squares optimization.
"""

from importlib.metadata import version
import math

import numpy as np
from scipy.signal import resample_poly


class ScalarAudioLoss:
    units = "native objective"

    def residual(self, samples, regions=range(5)):
        if tuple(regions) != tuple(range(5)):
            raise ValueError("Perceptual comparison uses the complete fixed event")
        return np.array([math.sqrt(max(0, self.score(samples)))])

    def diagnostics(self, samples):
        return dict(score=self.score(samples), units=self.units)


class AuralossMel(ScalarAudioLoss):
    """Library mel MR-STFT: spectral convergence plus log-magnitude L1."""

    def __init__(self, reference, rate, weighted=False):
        import torch
        from auraloss.freq import MultiResolutionSTFTLoss

        self.torch = torch
        self.target = torch.tensor(reference, dtype=torch.float32)[None, None]
        self.specification = dict(
            implementation="auraloss",
            version=version("auraloss"),
            fft_sizes=[512, 2048, 8192],
            hop_sizes=[128, 512, 2048],
            win_lengths=[512, 2048, 8192],
            scale="mel",
            n_bins=64,
            sample_rate=rate,
            perceptual_weighting=weighted,
            scale_invariance=False,
            w_sc=1.0,
            w_log_mag=1.0,
            eps=1e-10,
        )
        settings = {
            k: v
            for k, v in self.specification.items()
            if k not in ("implementation", "version")
        }
        self.loss = MultiResolutionSTFTLoss(**settings)

    def score(self, samples):
        if len(samples) != self.target.shape[-1]:
            raise ValueError("Candidate/reference lengths must match")
        with self.torch.inference_mode():
            value = self.torch.tensor(samples, dtype=self.torch.float32)[None, None]
            return float(self.loss(value, self.target))


class JtfsLoss(ScalarAudioLoss):
    """Log-compressed joint scattering, retaining time instead of global pooling.

    Uses the DAFx2022 authors' implementation. This is a direct feature distance,
    not a reproduction of PNP's trained inverse network or its parameter metric.
    """

    def __init__(self, reference, rate, device="cpu"):
        import torch
        from wavespin import TimeFrequencyScattering1D

        self.torch, self.rate, self.device = torch, rate, device
        self.length = len(reference)
        self.analysis_rate = 16000
        settings = dict(
            shape=round(len(reference) * self.analysis_rate / rate),
            J=10,
            Q=8,
            T=256,
            J_fr=3,
            Q_fr=1,
            F=4,
            average=True,
            average_fr=True,
            out_type="array",
            frontend="torch",
            pad_mode="zero",
        )
        self.transform = TimeFrequencyScattering1D(**settings)
        if device == "cuda":
            self.transform.cuda()
        raw = self.features(reference)
        self.floor = max(float(np.max(raw)) * 1e-3, 1e-10)
        self.target = np.log1p(raw / self.floor)
        self.specification = dict(
            implementation="wavespin",
            version=version("wavespin"),
            upstream_commit="5ff1b72785703cd09d4b1bf4b5f52ffcc8a926ae",
            settings=settings,
            analysis_rate=self.analysis_rate,
            compression="log1p(abs(coefficients)/reference_floor)",
            floor=self.floor,
            distance="mean squared feature difference",
            normalization="none",
            time_pooling="16ms; no global pooling",
        )

    def features(self, samples):
        if len(samples) != self.length:
            raise ValueError("Candidate/reference lengths must match")
        divisor = math.gcd(self.rate, self.analysis_rate)
        audio = resample_poly(
            samples, self.analysis_rate // divisor, self.rate // divisor
        )
        with self.torch.inference_mode():
            tensor = self.torch.tensor(
                audio, dtype=self.torch.float32, device=self.device
            )
            return np.abs(self.transform(tensor).cpu().numpy())

    def score(self, samples):
        return float(
            np.mean((np.log1p(self.features(samples) / self.floor) - self.target) ** 2)
        )
