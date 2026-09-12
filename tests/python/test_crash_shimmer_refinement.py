"""Shared crash edits preserve public controls, dormant handles and decay sparsity."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest

pytest.importorskip("torch")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from crash_refinement_search import (
    shared_edit,
    texture_cases,
    damping_cases,
    fit_upper_balance,
    tuning_decay_edit,
)
from triggerfish_percussion.saved_fit_renderer import SavedFitRenderer
from crash_texture_diagnostics import CrashBalance, ridge_contrast
from crash_refinement_contracts import training_seeds, audit_seeds, prepare_output


def parameters():
    fit = json.loads((ROOT / "workbench/web/crash_calibration.fit.json").read_text())
    return SavedFitRenderer.parameters(fit)


def test_shared_edit_identity_and_unchanged_control_surface():
    p = parameters()
    result = shared_edit(
        p,
        [
            1,
            0,
            p["body_decay_seconds_0"],
            p["body_decay_seconds_7"],
            p["bloom_rate"],
            p["body_brightness"],
            p["body_excitation_centre"],
            0,
            0,
        ],
    )
    assert result == p


def test_shared_pitch_scaling_and_shelves_do_not_touch_local_decay_or_gains():
    p = parameters()
    result = shared_edit(p, [1.05, 0, 20, 0.6, 5, -10, 2000, 2, -2])
    assert set(result) == set(p)
    for i in range(32):
        if p[f"resolved_level_{i}"] <= -71.99:
            assert result[f"resolved_frequency_{i}"] == p[f"resolved_frequency_{i}"]
            assert result[f"resolved_level_{i}"] == p[f"resolved_level_{i}"]
        else:
            assert result[f"resolved_frequency_{i}"] == pytest.approx(
                1.05 * p[f"resolved_frequency_{i}"]
            )
        assert result[f"resolved_turbulence_{i}"] == p[f"resolved_turbulence_{i}"]
    for key in ["direct_gain", "field_gain", "model_level_db", "body_excitation"]:
        assert result[key] == p[key]
    assert {k: v for k, v in result.items() if k.startswith("body_decay_active_")} == {
        k: v for k, v in p.items() if k.startswith("body_decay_active_")
    }


def test_texture_trials_change_no_pitch_decay_or_output_gains():
    p = parameters()
    for name, trial in texture_cases(p):
        assert set(trial) == set(p)
        for key in p:
            if key.startswith(("resolved_", "body_decay_")) or key in [
                "model_level_db",
                "field_gain",
                "direct_gain",
            ]:
                assert trial[key] == p[key], (name, key)
        if name != "unchanged":
            assert trial["output_eq_enabled"] == 0
        assert np.isfinite(list(trial.values())).all()


@pytest.mark.parametrize("values", [np.zeros(8), np.full(9, np.nan)])
def test_bad_shared_coordinates_fail(values):
    with pytest.raises(ValueError):
        shared_edit({}, values)


def test_texture_locked_edit_changes_only_shared_frequencies_and_endpoints():
    p = parameters()
    p["resolved_level_0"] = -71.5
    p["resolved_level_1"] = -71.98
    changed = tuning_decay_edit(p, [20, 1.4, 1.02, -0.01])
    for key in p:
        if key.startswith("resolved_frequency_") or key in (
            "body_decay_seconds_0",
            "body_decay_seconds_7",
        ):
            continue
        assert changed[key] == p[key], key
    assert changed["body_decay_seconds_0"] == 20
    assert changed["body_decay_seconds_7"] == 1.4


def test_upper_balance_is_a_broad_level_edit_only():
    p = parameters()
    p["resolved_level_0"] = -71.5
    p["resolved_level_20"] = -71.5
    trials = []
    fit_upper_balance(p, lambda values, name: trials.append(values))
    assert len(trials) == 6
    for trial in trials:
        for key in p:
            if not key.startswith("resolved_level_"):
                assert trial[key] == p[key]
        for i in range(32):
            key = f"resolved_level_{i}"
            if p[key] <= -71.99 or p[f"resolved_frequency_{i}"] <= 1500:
                assert trial[key] == p[key]
            else:
                assert p[key] - 6 <= trial[key] <= p[key]


def test_reference_self_fit_and_wrong_upper_decay_are_distinguished():
    rate = 16000
    t = np.arange(6 * rate) / rate
    noise = np.random.default_rng(91).normal(0, 0.1, t.size)
    body = 0.2 * np.sin(2 * np.pi * 125 * t) * np.exp(-t / 1.8)
    target = body + noise * np.exp(-t / 0.3)
    loss = CrashBalance(target, rate)
    identical = loss.components(target)
    slow = loss.components(body + noise * np.exp(-t / 0.9))
    quiet = loss.components(target * 0.25)
    assert identical["score"] < 1e-6
    assert slow["shape_error_db"] > 1
    assert slow["score"] > identical["score"] + 0.1
    assert quiet["score"] > identical["score"] + 0.1


def test_line_contrast_separates_stable_ridges_from_noise():
    rate = 16000
    t = np.arange(3 * rate) / rate
    noise = np.random.default_rng(7).normal(0, 1, len(t))
    ring = (
        sum(np.sin(2 * np.pi * f * t) for f in [1700, 1920, 2340, 2710]) + noise * 0.01
    )
    tonal = ridge_contrast(ring, rate)["flatness_db"][1][0]
    wash = ridge_contrast(noise, rate)["flatness_db"][1][0]
    assert tonal < wash - 5


def test_sparse_decay_changes_only_one_shared_knot_and_upper_endpoint():
    p = parameters()
    # This search explicitly requires two endpoints, independent of the live preset.
    p.update({key: 0 for key in p if key.startswith("body_decay_active_")})
    allowed = {
        "body_decay_active_1",
        "body_decay_frequency_1",
        "body_decay_seconds_1",
        "body_decay_seconds_7",
    }
    for _, trial in damping_cases(p):
        assert set(trial) == set(p)
        assert (
            sum(
                v >= 0.5 for k, v in trial.items() if k.startswith("body_decay_active_")
            )
            == 1
        )
        assert all(trial[k] == p[k] for k in p if k not in allowed)
    with pytest.raises(ValueError):
        list(damping_cases(dict(p, body_decay_active_1=1)))


def test_quiet_plateau_floor_does_not_discard_decaying_low_ring():
    rate = 16000
    t = np.arange(6 * rate) / rate
    rng = np.random.default_rng(341)
    body = 0.2 * np.sin(2 * np.pi * 125 * t) * np.exp(-t / 1.8)
    noise = rng.normal(0, 1, len(t)) * (0.1 * np.exp(-t / 0.2) + 0.00003)
    loss = CrashBalance(body + noise, rate)
    shape = loss.original.shape
    band = np.flatnonzero((shape.edges[:-1] < 125) & (shape.edges[1:] > 125))[0]
    assert loss.tail_floor[band] == -300
    assert np.any(loss.tail_floor > -300)
    assert not loss.late[:9].any()


def test_audit_seeds_never_overlap_training():
    for primary in (1982, 73519, 41273, 1396978464):
        train = training_seeds(primary)
        assert len(set(train)) == 2
        assert not set(train).intersection(audit_seeds(primary))
    assert audit_seeds(73519, 41273) == [69101, 91387]


def test_output_must_not_overwrite_previous_results(tmp_path):
    output = tmp_path / "new-fit"
    prepare_output(output)
    trace = output / "source.fit.json"
    trace.write_text("saved input", encoding="utf8")
    with pytest.raises(ValueError, match="fresh --output"):
        prepare_output(output)
    assert trace.read_text(encoding="utf8") == "saved input"


def test_locked_solver_start_is_inside_search_bounds(monkeypatch):
    import crash_refinement_search as search

    p = dict(parameters(), body_decay_seconds_0=0.5, body_decay_seconds_7=20)

    def optimizer(fun, start, **options):
        assert np.all((start >= 0) & (start <= 1))
        return fun(start)

    monkeypatch.setattr(search, "minimize", optimizer)
    search.fit_locked_texture(p, lambda values, stage: 0, 1)
    assert p["body_decay_seconds_0"] == 0.5
