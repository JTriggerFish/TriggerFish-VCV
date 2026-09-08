"""Band-envelope display must preserve candidate gain differences."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from triggerfish_percussion.audio_io import AudioBuffer, write_wav

spec = importlib.util.spec_from_file_location(
    "decay_plot", Path(__file__).parents[2] / "tools/plot_decay_comparison.py"
)
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)


def test_candidate_gain_is_not_normalized(tmp_path):
    rate = 48000
    time = np.arange(rate * 4) / rate
    reference = np.random.default_rng(71).normal(size=time.size) * np.exp(-3 * time)
    for name, samples in (("reference", reference), ("candidate", reference * 2)):
        write_wav(tmp_path / f"{name}.wav", AudioBuffer(samples, rate))
    result = plots.figure(tmp_path)
    assert len(result.data) == 12
    for index in range(0, 12, 2):
        ref, candidate = result.data[index : index + 2]
        selected = np.asarray(ref.x) < 1
        difference = np.asarray(candidate.y)[selected] - np.asarray(ref.y)[selected]
        assert np.allclose(difference, 20 * np.log10(2), atol=1e-5)


def test_mismatched_duration_is_rejected(tmp_path):
    for name, frames in (("reference", 48000 * 4), ("candidate", 48000 * 3)):
        write_wav(tmp_path / f"{name}.wav", AudioBuffer(np.zeros(frames), 48000))
    with pytest.raises(ValueError, match="durations differ"):
        plots.figure(tmp_path)
