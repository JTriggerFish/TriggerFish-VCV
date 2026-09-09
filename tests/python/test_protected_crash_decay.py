"""Decay analysis preserves a user's front; it never auto-normalizes playback."""

import sys
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip("torch")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from refine_user_crash_decay import ProtectedDecay


@pytest.fixture(scope="module")
def case():
    rate = 24000
    time = np.arange(6 * rate) / rate
    audio = np.random.default_rng(54).normal(0, 0.1, len(time)) * 10 ** (-3 * time / 8)
    return audio, ProtectedDecay(audio, audio, rate)


def test_identical_signal_has_zero_error(case):
    audio, loss = case
    assert np.linalg.norm(loss.residual(audio)) < 1e-10
    assert loss.diagnostics(audio)["front_rms_db"] < 1e-10


def test_gain_change_is_not_hidden_by_decay_shape_alignment(case):
    audio, loss = case
    diagnostics = loss.diagnostics(audio * 2)
    assert diagnostics["front_rms_db"] > 5.9
    assert diagnostics["score"] > 1


def test_tail_edit_is_measured_without_front_change(case):
    audio, loss = case
    changed = audio.copy()
    changed[2 * 24000 :] *= 0.3
    assert loss.diagnostics(changed)["front_rms_db"] < 1e-10
    assert np.linalg.norm(loss.residual(changed)) > 1
