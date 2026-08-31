"""Numerical analysis primitives for calibrated TriggerFish percussion."""

from .audio_io import AudioBuffer, read_wav, resample_audio, write_wav
from .erb import ErbFilterbank, erb_rate, frequency_from_erb_rate
from .transforms import StftConfig, StftResult, multiresolution_stft, stft

__all__ = [
    "AudioBuffer",
    "ErbFilterbank",
    "StftConfig",
    "StftResult",
    "erb_rate",
    "frequency_from_erb_rate",
    "multiresolution_stft",
    "read_wav",
    "resample_audio",
    "stft",
    "write_wav",
]
