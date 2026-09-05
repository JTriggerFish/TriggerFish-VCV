"""Reference proposals cover bass and upper ringing, without hidden mode rules."""

import numpy as np
from triggerfish_percussion.modal_fit_initialization import (
    spectral_mode_candidates,
    reference_modal_starts,
)


def test_proposals_reach_upper_ringing():
    rate = 48000
    time = np.arange(rate) / rate
    signal = np.exp(-8 * time) * (
        np.sin(2 * np.pi * 70 * time) + 0.3 * np.sin(2 * np.pi * 1200 * time)
    )
    peaks = spectral_mode_candidates(signal, rate)
    assert any(abs(item["frequency"] - 70) < 6 for item in peaks)
    assert any(abs(item["frequency"] - 1200) < 6 for item in peaks)
    starts = reference_modal_starts(
        {"equalizer_mode": 0, "resonance_decay_seconds": 0.6}, peaks
    )
    assert any(values["resonance_frequency_1"] > 1000 for _, values in starts)
    assert all(values["resonance_level_0"] == 0 for _, values in starts)
    assert all(values["resonance_decay_seconds"] == 0.6 for _, values in starts)
    assert all(values["equalizer_mode"] == 0 for _, values in starts)


def test_silence_creates_no_fake_modes():
    assert spectral_mode_candidates(np.zeros(48000), 48000) == []
    assert reference_modal_starts({}, []) == []
