"""Versioned on-disk cache for canonical STFT results."""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .audio_io import AudioBuffer
from .transforms import StftConfig, StftResult, TRANSFORM_VERSION, stft


class TransformCache:
    def __init__(self, directory: str | Path):
        self.directory = Path(directory)

    def stft(self, audio: AudioBuffer, config: StftConfig) -> StftResult:
        key = _cache_key(audio, config)
        path = self.directory / f"{key}.npz"
        if path.exists():
            return _load(path, config)
        result = stft(audio.mono().samples, audio.sample_rate, config)
        self.directory.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".tmp.npz")
        np.savez_compressed(
            temporary,
            frequencies_hz=result.frequencies_hz,
            times_seconds=result.times_seconds,
            spectrum=result.spectrum,
        )
        temporary.replace(path)
        return result


def _cache_key(audio: AudioBuffer, config: StftConfig) -> str:
    payload = {
        "audio": audio.content_hash(),
        "config": config.__dict__,
        "version": TRANSFORM_VERSION,
    }
    return sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _load(path: Path, config: StftConfig) -> StftResult:
    with np.load(path, allow_pickle=False) as stored:
        return StftResult(
            stored["frequencies_hz"],
            stored["times_seconds"],
            stored["spectrum"],
            config,
        )
