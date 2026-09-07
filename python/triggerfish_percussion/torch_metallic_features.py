"""Optional differentiable analysis of audio, not a percussion DSP surrogate."""

import numpy as np
from scipy.fft import next_fast_len
from scipy.signal import sosfilt
import torch
import torch.nn.functional as functional

from .power_envelope import _kernel
from .metallic_balance_features import ATTACK_BINS


class TorchMetallicFeatures:
    def __init__(self, features):
        self.spec = features
        self.frames, self.rate = features.frames, features.rate
        # A full second covers the slowest 40Hz band's IIR impulse to far below
        # float64 measurement relevance. Validation checks the actual features.
        impulse = np.zeros(round(self.rate))
        impulse[0] = 1
        kernels = np.array([sosfilt(sos, impulse) for sos in features.filters])
        self.filter_length = next_fast_len(self.frames + len(impulse) - 1)
        self.filters = torch.fft.rfft(torch.tensor(kernels), n=self.filter_length)
        self.radius, gaussian = _kernel(0.012 * self.rate)
        self.power_length = next_fast_len(self.frames + 4 * self.radius)
        self.gaussian = torch.fft.rfft(torch.tensor(gaussian), n=self.power_length)
        self.reflection = torch.tensor(
            np.pad(np.arange(self.frames), (self.radius, self.radius), mode="symmetric")
        )

    @staticmethod
    def stft(samples, window, hop):
        weights = torch.hann_window(window, periodic=True, dtype=samples.dtype)
        transformed = torch.stft(
            samples,
            window,
            hop,
            window=weights,
            center=True,
            pad_mode="constant",
            return_complex=True,
        )
        scale = torch.full((window // 2 + 1,), 2 / weights.sum(), dtype=samples.dtype)
        scale[0] *= 0.5
        scale[-1] *= 0.5
        return transformed * scale[:, None]

    @staticmethod
    def smooth_frequency(values, sigma):
        radius, weights = _kernel(sigma)
        indices = torch.tensor(
            np.pad(np.arange(values.shape[-1]), (radius, radius), mode="symmetric")
        )
        padded = values[:, indices].unsqueeze(1)
        return functional.conv1d(padded, torch.tensor(weights)[None, None, :]).squeeze(
            1
        )

    def __call__(self, samples):
        if (
            samples.device.type != "cpu"
            or samples.is_complex()
            or samples.shape != (self.frames,)
            or not torch.isfinite(samples).all()
        ):
            raise ValueError(
                "Expected finite real mono CPU audio of the reference duration"
            )
        # Analysis constants are float64. Casting here also preserves autograd
        # back to a float32 audio producer rather than failing in convolution.
        samples = samples.to(dtype=torch.float64)
        spectrum = torch.fft.rfft(samples, n=self.filter_length)
        filtered = torch.fft.irfft(self.filters * spectrum, n=self.filter_length)[
            :, : self.frames
        ]
        padded = filtered[:, self.reflection].square()
        power = torch.fft.irfft(
            torch.fft.rfft(padded, n=self.power_length) * self.gaussian,
            n=self.power_length,
        )
        envelope = power[
            :, 2 * self.radius : 2 * self.radius + self.frames : self.spec.hop
        ].clamp_min(0)
        transformed = self.stft(samples, 16384, 4096)
        frequencies = torch.arange(transformed.shape[0]) * self.rate / 16384
        selected = (frequencies >= 40) & (frequencies <= min(16000, 0.49 * self.rate))
        power = transformed[selected].abs().square()
        times = torch.arange(power.shape[1]) * 4096 / self.rate
        spectra = torch.stack(
            [
                power[:, (times >= a) & (times < b)].mean(dim=1)
                for a, b in self.spec.regions
            ]
        )
        spectra = self.smooth_frequency(spectra, 0.5)
        attack = self.stft(samples[: round(0.14 * self.rate)], 512, 128)
        frequencies = torch.arange(attack.shape[0]) * self.rate / 512
        selected = (frequencies >= 80) & (frequencies <= min(16000, 0.49 * self.rate))
        times = torch.arange(attack.shape[1]) * 128 / self.rate
        return dict(
            envelope=envelope,
            spectrum=spectra,
            attack=attack[selected][:, times < 0.12].abs().square(),
            transient=torch.stack(
                [
                    samples[round(a * self.rate) : round(b * self.rate)].square().mean()
                    for a, b in ATTACK_BINS
                ]
            ),
        )
