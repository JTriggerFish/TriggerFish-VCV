"""Document broad-band loss ambiguities and the independent texture diagnostic."""

import importlib
from pathlib import Path

import pytest
import numpy as np

pytest.importorskip("torch")
pytest.importorskip("auraloss")


def test_sizzle_diagnostic_exposes_broad_power_blind_spots(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[2] / "tools"))
    audit = importlib.import_module("audit_gong_sizzle")
    result = audit.blind_spots(44100)
    assert result["within_band_pitch_shift_8_to_13khz_broad_envelope_error_db"] < 0.01
    assert result["forty_hz_am_broad_envelope_error_db"] < 0.01
    assert result["forty_hz_am_texture_distance"] > 1


def test_ridge_diagnostic_distinguishes_noise_from_resolved_partials(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[2] / "tools"))
    audit = importlib.import_module("audit_gong_sizzle")
    rate = 32000
    time = np.arange(2 * rate) / rate
    tones = sum(np.sin(2 * np.pi * f * time) for f in (4000, 5000, 9000, 11000))
    noise = np.random.default_rng(7).normal(size=time.size)
    assert np.all(
        np.array(audit.ridge_contrast(tones, rate))
        > np.array(audit.ridge_contrast(noise, rate)) + 1
    )


def test_sizzle_colour_preserves_low_body_and_modal_dynamics(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[2] / "tools"))
    polish = importlib.import_module("polish_gong_sizzle")
    base = {"bloom_rate": 6, "model_level_db": 0, "body_decay_seconds_0": 10}
    for i, f in enumerate(np.geomspace(120, 15000, 32)):
        base[f"resolved_frequency_{i}"] = f
        base[f"resolved_level_{i}"] = -10
    base["resolved_level_20"] = -72
    before = dict(base)
    result = polish.colour(base, -6)
    assert base == before
    for key, value in base.items():
        if not key.startswith("resolved_level_"):
            assert result[key] == value
    for i in range(32):
        key = f"resolved_level_{i}"
        if (
            base[f"resolved_frequency_{i}"] <= 900
            or base[f"resolved_frequency_{i}"] >= 12000
        ):
            assert result[key] == base[key]
        assert result[key] <= base[key]
    assert result["resolved_level_20"] == -72
