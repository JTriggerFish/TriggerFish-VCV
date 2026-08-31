from pathlib import Path

import numpy as np
import pytest
from scipy.io import wavfile

from triggerfish_percussion.alignment import (
    detect_impact_onset,
    measure_alignment,
    shift_with_zeros,
)
from triggerfish_percussion.audio_io import (
    AudioBuffer,
    read_wav,
    resample_audio,
    write_wav,
)
from triggerfish_percussion.segmentation import nominal_regions


def test_float_wav_round_trip_preserves_level_and_channels(tmp_path: Path):
    samples = np.column_stack(
        (np.linspace(-0.8, 0.8, 257), np.linspace(0.4, -0.4, 257))
    )
    path = tmp_path / "level.wav"
    write_wav(path, AudioBuffer(samples, 48_000))
    restored = read_wav(path)
    assert restored.sample_rate == 48_000
    assert restored.samples == pytest.approx(samples, abs=3.0e-8)
    assert restored.mono("left").samples == pytest.approx(samples[:, 0])
    assert restored.mono("right").samples == pytest.approx(samples[:, 1])


def test_integer_wav_conversion_uses_full_scale(tmp_path: Path):
    path = tmp_path / "integer.wav"
    wavfile.write(path, 8_000, np.array([-32768, 0, 32767], dtype=np.int16))
    samples = read_wav(path).samples
    assert samples[0] == -1.0
    assert samples[1] == 0.0
    assert samples[2] == pytest.approx(32767 / 32768)


def test_resampling_preserves_duration_and_sine_level():
    time = np.arange(4_800) / 48_000
    source = AudioBuffer(0.5 * np.sin(2 * np.pi * 997 * time), 48_000)
    result = resample_audio(source, 44_100)
    assert result.samples.size == 4_410
    assert np.sqrt(np.mean(result.samples[100:-100] ** 2)) == pytest.approx(
        0.5 / np.sqrt(2), rel=2.0e-3
    )


def test_impact_onset_and_pair_alignment_are_sample_preserving():
    sample_rate = 48_000
    reference = np.zeros(8_000)
    onset = 2_000
    reference[onset : onset + 16] = np.hanning(32)[16:]
    candidate = shift_with_zeros(reference, 37)
    assert detect_impact_onset(reference, sample_rate) == pytest.approx(onset, abs=48)
    alignment = measure_alignment(reference, candidate, sample_rate)
    assert alignment.candidate_lag_samples == pytest.approx(-37, abs=2)
    aligned = shift_with_zeros(candidate, alignment.candidate_lag_samples)
    assert aligned == pytest.approx(reference)


def test_nominal_regions_are_named_and_nonoverlapping():
    regions = nominal_regions(100, 48_100, 48_000)
    slices = regions.slices()
    assert tuple(slices) == ("contact", "bloom", "early_body", "tail")
    assert slices["contact"].start == 100
    assert slices["tail"].stop == 48_100
