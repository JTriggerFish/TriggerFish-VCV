"""A cropped difference view cannot include a later impulse."""

import importlib.util
from pathlib import Path
import numpy as np
from triggerfish_percussion.audio_io import AudioBuffer, write_wav

spec = importlib.util.spec_from_file_location(
    "difference_plot", Path(__file__).parents[2] / "tools/plot_spectral_difference.py"
)
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)


def test_crop_precedes_transform(tmp_path):
    rate = 48000
    reference = np.zeros(rate)
    candidate = reference.copy()
    candidate[round(0.125 * rate)] = 1
    for name, samples in (("reference", reference), ("candidate", candidate)):
        write_wav(tmp_path / f"{name}.wav", AudioBuffer(samples, rate))
    cropped = plots.figure(tmp_path, window=512, hop=64, seconds=0.12)
    full = plots.figure(tmp_path, window=512, hop=64)
    assert not np.asarray(cropped.data[2].z).any()
    assert np.asarray(full.data[2].z).any()
