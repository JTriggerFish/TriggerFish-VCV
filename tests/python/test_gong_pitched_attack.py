"""Attack studies do not quietly retune the upper series, damping or gains."""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("torch")
pytest.importorskip("auraloss")


@pytest.mark.parametrize("study", ["texture", "low-core", "five-core", "front"])
def test_attack_studies_freeze_unrelated_parameters(monkeypatch, study):
    root = Path(__file__).parents[2]
    monkeypatch.syspath_prepend(str(root / "tools"))
    module = importlib.import_module("refine_gong_pitched_attack")
    fit = json.loads(
        (root / "workbench/web/gong_calibration.fit.json").read_text(encoding="utf8")
    )
    base = {
        k: v
        for node in fit["instrument"]["nodes"]
        for k, v in node["parameters"].items()
    }
    original = dict(base)
    rows = list(module.candidates(base, study))
    assert len({name for name, _ in rows}) == len(rows)
    for _, p in rows:
        assert p.keys() == base.keys()
        assert np.isfinite(list(p.values())).all()
        for k in base:
            permitted = k in ("bloom_rate", "bloom_energy_acceleration") or any(
                k == f"resolved_{field}_{i}"
                for field in ("frequency", "level", "turbulence")
                for i in range(5)
            )
            if not permitted:
                assert p[k] == base[k]
        if study == "texture":
            for i in range(32):
                assert p[f"resolved_frequency_{i}"] == base[f"resolved_frequency_{i}"]
    assert base == original


def test_attack_guard_handles_exact_baseline_without_normalizing(monkeypatch):
    pytest.importorskip("auraloss")
    import torch

    torch.set_num_threads(1)
    root = Path(__file__).parents[2]
    monkeypatch.syspath_prepend(str(root / "tools"))
    module = importlib.import_module("polish_gong_pitched_attack")
    rate = 24000
    time = np.arange(6 * rate) / rate
    ref = 0.1 * np.sin(2 * np.pi * 344 * time) * np.exp(-time)
    loss = module.PitchedAttackLoss(ref, rate, [ref])
    assert loss.score(ref) == 0
    assert np.isfinite(loss.score(ref * 2)) and loss.score(ref * 2) > 0
