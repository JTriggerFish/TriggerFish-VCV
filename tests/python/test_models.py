import numpy as np
import pytest

import _triggerfish_dsp as dsp


def test_vca_binding_is_finite_and_length_preserving():
    audio = np.sin(np.linspace(0.0, 8.0 * np.pi, 1024, dtype=np.float32))
    cv = np.linspace(0.0, 1.0, audio.size, dtype=np.float32)
    output = dsp.vca_transistor(audio, cv, 48_000.0)
    assert output.shape == audio.shape
    assert np.isfinite(output).all()


def test_vdpo_binding_is_finite_and_length_preserving():
    size = 256
    audio = np.zeros(size)
    damping = np.full(size, 0.5)
    angular_frequency = np.full(size, 2.0 * np.pi * 261.625565)
    output = dsp.vdpo(audio, damping, angular_frequency, 48_000.0)
    assert output.shape == audio.shape
    assert np.isfinite(output).all()


def test_bindings_reject_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        dsp.vca_transistor(np.zeros(4), np.zeros(3), 48_000.0)
