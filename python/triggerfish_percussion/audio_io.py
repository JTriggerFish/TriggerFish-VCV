"""Level-preserving WAV I/O and deterministic resampling."""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly


@dataclass(frozen=True)
class AudioBuffer:
    samples: np.ndarray
    sample_rate: int

    def __post_init__(self) -> None:
        samples = np.asarray(self.samples, dtype=np.float64)
        if samples.ndim not in (1, 2) or samples.size == 0:
            raise ValueError("audio must be a non-empty mono or frame-by-channel array")
        if self.sample_rate < 1 or not np.isfinite(samples).all():
            raise ValueError("audio sample rate and samples must be finite")
        object.__setattr__(self, "samples", samples)

    def mono(self, policy: str = "mean") -> "AudioBuffer":
        if self.samples.ndim == 1:
            return self
        if policy == "mean":
            samples = np.mean(self.samples, axis=1)
        elif policy == "left":
            samples = self.samples[:, 0]
        elif policy == "right":
            samples = self.samples[:, -1]
        else:
            raise ValueError(f"unknown channel policy: {policy}")
        return AudioBuffer(samples, self.sample_rate)

    def content_hash(self) -> str:
        digest = sha256()
        digest.update(str(self.sample_rate).encode("ascii"))
        digest.update(np.ascontiguousarray(self.samples).view(np.uint8))
        return digest.hexdigest()


def _float_samples(samples: np.ndarray) -> np.ndarray:
    if np.issubdtype(samples.dtype, np.floating):
        return samples.astype(np.float64)
    if np.issubdtype(samples.dtype, np.signedinteger):
        scale = float(
            max(abs(np.iinfo(samples.dtype).min), np.iinfo(samples.dtype).max)
        )
        return samples.astype(np.float64) / scale
    if np.issubdtype(samples.dtype, np.unsignedinteger):
        info = np.iinfo(samples.dtype)
        midpoint = (info.max + 1) / 2
        return (samples.astype(np.float64) - midpoint) / midpoint
    raise TypeError(f"unsupported WAV dtype: {samples.dtype}")


def read_wav(path: str | Path, channel_policy: str | None = None) -> AudioBuffer:
    sample_rate, samples = wavfile.read(Path(path), mmap=False)
    audio = AudioBuffer(_float_samples(samples), int(sample_rate))
    return audio if channel_policy is None else audio.mono(channel_policy)


def write_wav(path: str | Path, audio: AudioBuffer) -> None:
    samples = np.asarray(audio.samples, dtype=np.float32)
    wavfile.write(Path(path), audio.sample_rate, samples)


def resample_audio(audio: AudioBuffer, sample_rate: int) -> AudioBuffer:
    if sample_rate < 1:
        raise ValueError("target sample rate must be positive")
    if sample_rate == audio.sample_rate:
        return audio
    divisor = np.gcd(audio.sample_rate, sample_rate)
    samples = resample_poly(
        audio.samples,
        sample_rate // divisor,
        audio.sample_rate // divisor,
        axis=0,
        padtype="line",
    )
    return AudioBuffer(samples, sample_rate)
