"""Gong experiments retain a complete visible surface and frozen controls."""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("torch")


@pytest.fixture
def trials(monkeypatch):
    root = Path(__file__).parents[2]
    monkeypatch.syspath_prepend(str(root / "tools"))
    module = importlib.import_module("refine_gong_slow_beating")
    fit = json.loads(
        (root / "workbench/web/gong_calibration.fit.json").read_text(encoding="utf8")
    )
    base = {
        key: value
        for node in fit["instrument"]["nodes"]
        for key, value in node["parameters"].items()
    }
    return module, base


@pytest.mark.parametrize(
    "study", ["ring-balance", "ring-texture", "low-core", "clean-ring"]
)
def test_trials_preserve_surface_gain_damping_and_upper_series(trials, study):
    module, base = trials
    original = dict(base)
    rows = list(module.variants(base, study))
    assert len({name for name, _ in rows}) == len(rows)
    for _, candidate in rows:
        assert set(candidate) == set(base)
        assert np.isfinite(list(candidate.values())).all()
        for key in base:
            if key.startswith(("output_", "body_decay_", "impact_")) or key in (
                "model_level_db",
                "body_excitation",
                "field_gain",
                "direct_gain",
            ):
                assert candidate[key] == base[key]
        for i in range(4, 32):
            for field in ("frequency", "level", "turbulence", "allocation"):
                key = f"resolved_{field}_{i}"
                if key in base:
                    assert candidate[key] == base[key]
    assert base == original


def test_shared_loss_prefers_exact_signal_and_retains_gain(trials):
    pytest.importorskip("auraloss")
    import torch
    from refine_gong_shared_envelope import GongEnvelopeLoss

    torch.set_num_threads(1)
    rate = 16000
    time = np.arange(rate * 6) / rate
    audio = (
        0.2 * np.sin(2 * np.pi * 120 * time) + 0.1 * np.sin(2 * np.pi * 540 * time)
    ) * np.exp(-time)
    loss = GongEnvelopeLoss(audio, rate)
    assert loss.score(audio) < 1e-6
    assert loss.score(audio * 1.4) > 0.1
