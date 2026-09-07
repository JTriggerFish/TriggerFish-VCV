"""Opt-in integration tests for published objective implementations."""

import numpy as np
import pytest

pytest.importorskip("torch")
pytest.importorskip("auraloss")
from triggerfish_percussion.perceptual_fit_losses import AuralossMel, JtfsLoss
from triggerfish_percussion.loss_perturbations import fault_examples


def signal(rate=16000, frequency=65):
    t = np.arange(round(1.2 * rate)) / rate
    return 0.2 * np.sin(2 * np.pi * frequency * t) * np.exp(-12 * t)


@pytest.mark.parametrize("weighted", [False, True])
def test_mel_identity_gain_and_added_tail(weighted):
    import torch

    torch.set_num_threads(2)
    audio = signal()
    loss = AuralossMel(audio, 16000, weighted)
    assert loss.score(audio) == 0
    assert loss.score(audio * 0.5) > 0.01
    # Isolated pitch change: identical duration, envelope and amplitude.
    assert loss.score(signal(frequency=75)) > loss.score(signal(frequency=67)) > 0
    noisy = audio.copy()
    noisy[8000:9000] += 0.03 * np.random.default_rng(4).standard_normal(1000)
    assert loss.score(noisy) > 0.01
    with pytest.raises(ValueError):
        loss.score(audio[:-1])


def test_jtfs_identity_gain_and_shape():
    pytest.importorskip("wavespin")
    audio = signal()
    loss = JtfsLoss(audio, 16000)
    assert loss.score(audio) == 0
    assert loss.score(audio * 0.5) > 0.001
    assert loss.score(signal(frequency=75)) > loss.score(signal(frequency=67)) > 0
    with pytest.raises(ValueError):
        loss.score(audio[:-1])


def test_perturbations_have_consistent_lengths_and_severity():
    audio = signal()
    rows = list(fault_examples(audio, 16000))
    assert len(rows) == 12
    for name, severity, value in rows:
        assert value.shape == audio.shape
        assert np.isfinite(value).all()
        assert not np.array_equal(value, audio)


def test_missing_bass_never_boosts_spectral_bins():
    audio = np.random.default_rng(8).standard_normal(19200)
    original = np.abs(np.fft.rfft(audio))
    previous = original
    for name, severity, value in fault_examples(audio, 16000):
        if name != "missing-bass":
            continue
        actual = np.abs(np.fft.rfft(value))
        assert np.all(actual <= previous + 1e-10)
        assert actual[:100].sum() < previous[:100].sum()
        previous = actual
