"""Canonical time-frequency transforms shared by losses and diagnostics."""

from dataclasses import dataclass

import numpy as np
from scipy.signal import ShortTimeFFT, windows

TRANSFORM_VERSION = "stft-v1"


@dataclass(frozen=True)
class StftConfig:
    window_samples: int
    hop_samples: int
    fft_samples: int | None = None

    def validate(self) -> None:
        fft_samples = self.fft_samples or self.window_samples
        if self.window_samples < 4 or self.hop_samples < 1:
            raise ValueError("STFT window and hop must be positive")
        if fft_samples < self.window_samples:
            raise ValueError("STFT FFT size must cover its window")


@dataclass(frozen=True)
class StftResult:
    frequencies_hz: np.ndarray
    times_seconds: np.ndarray
    spectrum: np.ndarray
    config: StftConfig

    @property
    def magnitude(self) -> np.ndarray:
        return np.abs(self.spectrum)

    @property
    def power(self) -> np.ndarray:
        return np.square(self.magnitude)

    def magnitude_db(self, floor_db: float = -120.0) -> np.ndarray:
        peak = max(float(np.max(self.magnitude)), np.finfo(float).tiny)
        floor = peak * 10.0 ** (floor_db / 20.0)
        return 20.0 * np.log10(np.maximum(self.magnitude, floor) / peak)


def stft(samples: np.ndarray, sample_rate: float, config: StftConfig) -> StftResult:
    config.validate()
    signal = np.asarray(samples, dtype=np.float64)
    if signal.ndim != 1 or signal.size == 0 or not np.isfinite(signal).all():
        raise ValueError("STFT input must be finite non-empty mono audio")
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise ValueError("STFT sample rate must be positive")
    window = windows.hann(config.window_samples, sym=False)
    transform = ShortTimeFFT(
        window,
        config.hop_samples,
        sample_rate,
        fft_mode="onesided2X",
        mfft=config.fft_samples or config.window_samples,
        scale_to="magnitude",
    )
    spectrum = transform.stft(signal)
    times = transform.t(signal.size)
    selected = (times >= 0.0) & (times <= signal.size / sample_rate)
    return StftResult(
        transform.f.copy(), times[selected], spectrum[:, selected], config
    )


def multiresolution_stft(
    samples: np.ndarray, sample_rate: float, configs: tuple[StftConfig, ...]
) -> tuple[StftResult, ...]:
    if not configs:
        raise ValueError("at least one STFT resolution is required")
    return tuple(stft(samples, sample_rate, config) for config in configs)
