import numpy as np

from triggerfish_percussion.reference_onset_audit import audit_onset


def test_quiet_preroll_flagged_without_automatic_alignment():
    rate = 48000
    x = np.random.default_rng(42).normal(size=rate // 4) * 0.0001
    t = np.arange(round(0.05 * rate)) / rate
    start = round(0.064 * rate)
    x[start : start + len(t)] += 0.3 * np.cos(2 * np.pi * 1100 * t) * np.exp(-t / 0.02)
    audit = audit_onset(x, rate, 0.05)
    assert audit["possible_flat_preroll"]
    assert abs(audit["clear_rise_offset_seconds"] - 0.014) < 1 / rate
    assert not audit["automatic_alignment"]
    assert not audit_onset(x, rate, 0.064)["possible_flat_preroll"]


def test_trimmed_gradual_contact_not_mislabelled_as_preroll():
    rate = 48000
    t = np.arange(rate // 4) / rate
    x = t * np.exp(-t / 0.01) * np.cos(2 * np.pi * 500 * t)
    assert not audit_onset(x, rate, 0)["possible_flat_preroll"]
