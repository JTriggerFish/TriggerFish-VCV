from pathlib import Path
import sys

import numpy as np
import pytest

import _triggerfish_dsp as native
from triggerfish_percussion.audio_io import AudioBuffer, read_wav
from triggerfish_percussion.crash_corpus import (
    isolate_capture,
    load_fit_cells,
    write_cells_manifest,
)
from triggerfish_percussion.crash_fitting import (
    CrashFitCell,
    _maximum_late_regrowth_db,
    causal_audio_quality,
    causal_feature_loss,
    causal_prefix_features,
    causal_prefix_quality,
    parameter_influence_diagnostics,
    perceptual_features,
    prefix_gate_passes,
    reference_modal_residual,
)
from triggerfish_percussion.crash_fit_causal import (
    FIRST_100MS_TRADEOFF_POLICY,
    STRICT_CAUSAL_POLICY,
    _feasible_local_population,
    _initial_population,
    _protected_quality_passes,
    _validated_parameter_start,
    _validate_frozen_stages,
    fit_causal_model,
)
from triggerfish_percussion.crash_fit_spectral_profile import (
    temporal_spectral_parameter_names,
)
from triggerfish_percussion.crash_fit_parameters import (
    ATTACK_PARAMETERS,
    CausalStage,
    FitParameter,
    SCREENED_ATTACK_STAGES,
    SCREENED_INITIAL_DECAY_STAGES,
    SINGLE_HIT_ATTACK_PARAMETERS,
    SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS,
    UNIFIED_CAUSAL_STAGES,
    fit_parameter_value,
    replace_fit_parameters,
    single_hit_stages,
)
from triggerfish_percussion.crash_fit_texture import (
    compare_texture_descriptors,
    texture_descriptor,
)
from triggerfish_percussion.crash_model import (
    CrashEvent,
    CrashFit,
    render_crash,
    render_crash_components,
    render_crash_sequence,
)
from triggerfish_percussion.report import (
    ReportPair,
    _common_magnitude_db,
    velocity_grid_audition_gain,
    write_comparison_report,
)

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
import build_private_cymbal_sweep as sweep  # noqa: E402


def test_native_crash_parameters_round_trip_and_render():
    fit = CrashFit(
        contact_noise_gain=1.7,
        contact_chirp_frequency_scale=1.2,
        strength_gamma=2.0,
        output_gain=2.0,
    )
    first = render_crash(fit, 0.1, strength=0.7, seed=12)
    repeated = render_crash(fit, 0.1, strength=0.7, seed=12)
    assert first.shape == (4800,)
    assert np.array_equal(first, repeated)
    assert np.isfinite(first).all()
    assert np.max(np.abs(first)) > 0


def test_python_modal_grid_matches_the_native_constructive_defaults():
    fit = CrashFit()
    parameters = native.CrashCymbalFitParameters()
    for field in (
        "sparse_frequency_hz",
        "sparse_decay_ratio",
        "sparse_amplitude",
        "body_decay_frequency_hz",
        "body_decay_seconds",
        "body_decay_active",
    ):
        np.testing.assert_allclose(getattr(fit, field), getattr(parameters, field))


def test_legacy_decay_upgrade_preserves_all_five_interior_knots():
    frequencies = (150.0, 500.0, 1500.0, 6000.0, 16000.0)
    seconds = (3.5, 2.8, 2.5, 4.5, 0.35)
    fit = CrashFit(
        body_decay_frequency_hz=frequencies,
        body_decay_seconds=seconds,
    )
    assert fit.body_decay_frequency_hz[1:6] == frequencies
    assert fit.body_decay_seconds[1:6] == seconds
    assert fit.body_decay_seconds[0] == seconds[0]
    assert fit.body_decay_seconds[-1] == seconds[-1]
    assert fit.body_decay_active == (True, True, True, True, True, True, False, True)


def test_dense_mode_seed_is_repeatable_and_object_specific():
    first = render_crash(CrashFit(dense_mode_seed=123), 0.1, seed=12)
    repeated = render_crash(CrashFit(dense_mode_seed=123), 0.1, seed=12)
    different = render_crash(CrashFit(dense_mode_seed=124), 0.1, seed=12)
    assert np.array_equal(first, repeated)
    assert not np.array_equal(first, different)


def test_crash_binding_exposes_implement_and_brush_contact_spread():
    fit = CrashFit()
    brush_tap = render_crash(fit, 0.3, implement=0.0, contact_spread=0.0, seed=19)
    brush_sweep = render_crash(fit, 0.3, implement=0.0, contact_spread=1.0, seed=19)
    stick = render_crash(fit, 0.3, implement=1.0, contact_spread=1.0, seed=19)
    assert not np.array_equal(brush_tap, brush_sweep)
    assert not np.array_equal(brush_sweep, stick)
    late = slice(round(0.1 * 48_000), round(0.25 * 48_000))
    assert np.linalg.norm(brush_sweep[late]) > np.linalg.norm(brush_tap[late])


def test_indexed_fit_parameter_updates_one_dense_envelope_node():
    initial = CrashFit()
    fit = replace_fit_parameters(initial, {"dense_gain_envelope_db[3]": 6.0})
    assert fit.dense_gain_envelope_db[3] == 6.0
    changed = np.flatnonzero(
        np.asarray(fit.dense_gain_envelope_db)
        != np.asarray(initial.dense_gain_envelope_db)
    )
    assert changed.tolist() == [3]
    assert fit_parameter_value(fit, "dense_gain_envelope_db[3]") == 6.0


def test_old_dense_gain_envelope_is_interpolated_to_current_grid():
    fit = CrashFit(dense_gain_envelope_db=(0.0, 2.0, 4.0, 6.0, 8.0, 10.0))
    assert len(fit.dense_gain_envelope_db) == 33
    assert fit.dense_gain_envelope_db[0] == 0.0
    assert fit.dense_gain_envelope_db[-1] == 10.0


def test_spectral_profile_refinement_cannot_collapse_source_balance():
    names = temporal_spectral_parameter_names()
    assert "direct_gain" not in names
    assert "sparse_gain" not in names
    assert "dense_gain" not in names


def test_unified_fit_can_identify_bloom_routing_separately_from_body_level():
    stages = {stage.name: stage for stage in UNIFIED_CAUSAL_STAGES}
    for name in ("unified-initial-body", "unified-bloom"):
        parameters = {item.name for item in stages[name].parameters}
        assert "bloom_body_gain" in parameters
        assert "field_gain" in parameters or name == "unified-bloom"


def test_screened_initial_decay_keeps_earlier_gates_and_one_final_gate():
    stages = SCREENED_INITIAL_DECAY_STAGES
    assert [stage.name for stage in stages[:3]] == [
        "impact-contact",
        "impact-balance",
        "contact-tail",
    ]
    assert not any(stage.requires_acceptance_gate for stage in stages[3:-1])
    assert stages[-1].requires_acceptance_gate
    names = [parameter.name for parameter in stages[-1].parameters]
    assert len(names) == len(set(names))
    for index in (1, 2, 3, 7):
        assert f"body_decay_seconds[{index}]" in names


def test_resume_skips_diagnostic_only_intermediate_quality_gate():
    diagnostic = CausalStage("diagnostic", 0.1, (), requires_acceptance_gate=False)
    _validate_frozen_stages((), CrashFit(), (diagnostic,), {}, {}, "next")


def test_intermediate_stage_must_preserve_prior_absolute_quality():
    stage = CausalStage("next", 0.1, (), requires_acceptance_gate=False)
    passing = causal_prefix_quality(
        causal_prefix_features(AudioBuffer(np.ones(720), 48_000), 0.015),
        causal_prefix_features(AudioBuffer(np.ones(720), 48_000), 0.015),
    )
    failing = type(passing)(
        envelope_rmse_db=passing.envelope_rmse_db,
        spectral_rmse_db=passing.spectral_rmse_db,
        spectral_p95_absolute_db=passing.spectral_p95_absolute_db,
        fine_spectrum_rmse_db=stage.maximum_fine_spectrum_rmse_db + 0.1,
    )
    assert _protected_quality_passes({(0, 0.015): passing}, stage, (0.015,))
    assert not _protected_quality_passes({(0, 0.015): failing}, stage, (0.015,))


def test_causal_population_scales_with_dimension_and_is_deterministic():
    start = np.full(18, 0.5)
    bounds = np.tile((0.0, 1.0), (18, 1))
    first = _initial_population(start, bounds, 1600, 42)
    second = _initial_population(start, bounds, 1600, 42)
    assert first.shape == (128, 18)
    assert np.array_equal(first, second)
    assert np.array_equal(first[0], start)
    assert np.all((first >= 0.0) & (first <= 1.0))


def test_causal_optimizer_rejects_out_of_bound_seed_instead_of_clipping():
    with pytest.raises(ValueError, match="outside"):
        _validated_parameter_start(
            np.asarray([2.0]), np.asarray([[0.0, 1.0]]), ("amount",)
        )


def test_feasible_population_stays_inside_a_hard_gate():
    start = np.zeros(4)
    bounds = np.tile((-1.0, 1.0), (4, 1))

    def objective(values):
        return 0.0 if np.linalg.norm(values) <= 0.15 else 1.0e6

    first, first_probes = _feasible_local_population(
        start, bounds, 16, objective, 1.0e6, 42
    )
    second, second_probes = _feasible_local_population(
        start, bounds, 16, objective, 1.0e6, 42
    )
    assert first.shape == (16, 4)
    assert np.array_equal(first, second)
    assert first_probes == second_probes
    assert np.all(np.linalg.norm(first, axis=1) <= 0.15)


def test_causal_prefix_features_detect_wrong_attack_spectrum():
    sample_rate = 48_000
    time = np.arange(round(0.015 * sample_rate)) / sample_rate
    envelope = np.exp(-300.0 * time)
    bright = AudioBuffer(envelope * np.sin(2.0 * np.pi * 9000.0 * time), sample_rate)
    dark = AudioBuffer(envelope * np.sin(2.0 * np.pi * 1200.0 * time), sample_rate)
    target = causal_prefix_features(bright, 0.015)
    assert causal_feature_loss(target, target) == 0.0
    assert causal_feature_loss(causal_prefix_features(dark, 0.015), target) > 1.0e-6


def test_silent_optimizer_candidate_has_finite_texture_descriptor():
    descriptor = texture_descriptor(AudioBuffer(np.zeros(720), 48_000), 0.015, 0)
    arrays = (
        descriptor.fine_spectrum_db,
        descriptor.centroid_log2_hz,
        descriptor.rolloff_log2_hz,
        descriptor.flatness_db,
        descriptor.crest_db,
    )
    assert all(array.size > 0 and np.isfinite(array).all() for array in arrays)
    assert np.isfinite(descriptor.ridge_ratio)


def test_texture_trajectories_keep_fixed_length_for_quiet_candidates():
    reference = AudioBuffer(np.random.default_rng(9).normal(size=4800), 48_000)
    silence = AudioBuffer(np.zeros(4800), 48_000)
    quality = compare_texture_descriptors(
        texture_descriptor(silence, 0.100, 0),
        texture_descriptor(reference, 0.100, 0),
    )
    assert all(np.isfinite(value) for value in quality.__dict__.values())


def test_causal_prefix_quality_reports_absolute_level_error():
    sample_rate = 48_000
    signal = np.random.default_rng(31).normal(size=720)
    target = causal_prefix_features(AudioBuffer(signal, sample_rate), 0.015)
    quieter = causal_prefix_features(AudioBuffer(0.1 * signal, sample_rate), 0.015)
    quality = causal_prefix_quality(quieter, target)
    assert quality.envelope_rmse_db == pytest.approx(20.0, abs=0.1)
    assert quality.envelope_maximum_absolute_db == pytest.approx(20.0, abs=0.1)
    assert quality.spectral_rmse_db == pytest.approx(20.0, abs=0.1)


def test_causal_prefix_gate_rejects_earlier_regression():
    baseline = {0.004: 1.0, 0.015: 3.0}
    improved = {0.004: 1.01, 0.015: 2.0}
    regressed = {0.004: 1.03, 0.015: 1.0}
    assert prefix_gate_passes(baseline, improved, (0.004,))
    assert not prefix_gate_passes(baseline, regressed, (0.004,))


def test_first_100ms_policy_counts_nested_attack_features_once():
    assert STRICT_CAUSAL_POLICY.weight(0.004) == 1.0
    assert FIRST_100MS_TRADEOFF_POLICY.weight(0.004) == 0.05
    assert FIRST_100MS_TRADEOFF_POLICY.weight(0.015) == 0.15
    assert FIRST_100MS_TRADEOFF_POLICY.weight(0.100) == 1.0
    assert FIRST_100MS_TRADEOFF_POLICY.protected_objective_tolerance == 0.35
    assert FIRST_100MS_TRADEOFF_POLICY.protected_acceptance_tolerance == 0.35
    assert FIRST_100MS_TRADEOFF_POLICY.require_absolute_gate
    assert FIRST_100MS_TRADEOFF_POLICY.stop_on_absolute_failure


def test_attack_audit_domain_exposes_contact_and_body_controls():
    names = {parameter.name for parameter in ATTACK_PARAMETERS}
    assert {"contact_noise_gain", "contact_chirp_gain", "dispersion_drive"} <= names
    assert {"dense_gain", "sparse_gain", "body_decay_seconds[0]"} <= names


def test_screened_attack_uses_soft_blocks_then_one_hard_joint_gate():
    assert len(SCREENED_ATTACK_STAGES) == 6
    assert all(
        not stage.requires_acceptance_gate for stage in SCREENED_ATTACK_STAGES[:-1]
    )
    assert SCREENED_ATTACK_STAGES[-1].requires_acceptance_gate
    joint_names = {item.name for item in SCREENED_ATTACK_STAGES[-1].parameters}
    assert "contact_chirp_gain" not in joint_names
    assert {"contact_noise_gain", "dense_gain", "dispersion_feedback"} <= joint_names


def test_single_hit_schedule_excludes_velocity_curve_controls():
    stages = single_hit_stages(SCREENED_ATTACK_STAGES)
    names = {parameter.name for stage in stages for parameter in stage.parameters}
    assert names.isdisjoint(SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS)
    assert {"contact_noise_gain", "dense_gain", "dispersion_feedback"} <= names
    audit_names = {parameter.name for parameter in SINGLE_HIT_ATTACK_PARAMETERS}
    assert audit_names.isdisjoint(SINGLE_HIT_UNIDENTIFIABLE_PARAMETERS)


def test_absolute_quality_rejects_tonal_substitute_for_noise_attack():
    sample_rate = 48_000
    sample_count = round(0.100 * sample_rate)
    time = np.arange(sample_count) / sample_rate
    envelope = np.exp(-18.0 * time)
    random = np.random.default_rng(91)
    reference = envelope * random.normal(size=sample_count)
    frequencies = np.geomspace(350.0, 15_000.0, 28)
    phases = random.uniform(0.0, 2.0 * np.pi, frequencies.size)
    tonal = np.sum(
        np.sin(2.0 * np.pi * frequencies[:, None] * time + phases[:, None]), axis=0
    )
    tonal *= envelope * np.sqrt(
        np.mean(reference**2) / np.mean((envelope * tonal) ** 2)
    )

    quality = causal_audio_quality(
        AudioBuffer(tonal, sample_rate),
        AudioBuffer(reference, sample_rate),
        0.100,
    )

    assert quality.flatness_rmse_db > 2.5
    assert quality.fine_spectrum_rmse_db > 4.0
    assert quality.ridge_ratio_absolute_error > 0.08


def test_parameter_influence_reports_contact_in_first_prefix():
    sample_rate = 48_000
    reference = AudioBuffer(np.zeros(sample_rate // 20), sample_rate)
    cell = CrashFitCell("influence", reference, 0.8, seed=17)
    diagnostics = parameter_influence_diagnostics(
        cell,
        CrashFit(),
        (FitParameter("contact_pulse_gain", 0.0, 3.0),),
        0.015,
    )
    parameter_result = diagnostics["parameters"][0]
    assert parameter_result["earliest_detectable_seconds"] == 0.0


def test_causal_fit_recovers_a_synthetic_attack_parameter():
    target = CrashFit(contact_chirp_frequency_scale=1.45)
    cell = CrashFitCell(
        "synthetic",
        AudioBuffer(render_crash(target, 0.06, strength=0.8, seed=23), 48_000),
        0.8,
        seed=23,
    )
    stage = CausalStage(
        "impact",
        0.004,
        (FitParameter("contact_chirp_frequency_scale", 0.5, 2.0),),
    )
    fitted, diagnostics = fit_causal_model(
        (cell,),
        CrashFit(contact_chirp_frequency_scale=0.75),
        40,
        17,
        (stage,),
    )
    assert diagnostics["stages"][0]["accepted"]
    assert diagnostics["stages"][0]["composite_objective_improved"]
    assert abs(fitted.contact_chirp_frequency_scale - 1.45) < 0.1


def test_causal_fit_stops_when_active_prefix_is_not_matched():
    target = CrashFit(contact_chirp_frequency_scale=1.7)
    cell = CrashFitCell(
        "unreachable",
        AudioBuffer(render_crash(target, 0.06, strength=0.8, seed=23), 48_000),
        0.8,
        seed=23,
    )
    stages = (
        CausalStage(
            "wrong-control",
            0.004,
            (FitParameter("output_gain", 0.1, 4.0),),
            maximum_envelope_rmse_db=0.01,
            maximum_spectral_rmse_db=0.01,
            maximum_spectral_p95_db=0.01,
        ),
        CausalStage(
            "must-not-run",
            0.015,
            (FitParameter("contact_chirp_frequency_scale", 0.5, 2.0),),
        ),
    )
    _, diagnostics = fit_causal_model((cell,), CrashFit(), 24, 7, stages)
    assert not diagnostics["completed"]
    assert diagnostics["blocked_stage"] == "wrong-control"
    assert {stage["stage"] for stage in diagnostics["stages"]} == {"wrong-control"}
    stage = diagnostics["stages"][0]
    assert (
        stage["candidate_composite_objective"] <= stage["baseline_composite_objective"]
    )
    assert stage["candidate_composite_objective"] == min(
        stage["global_objective"], stage["local_objective"]
    )
    assert (
        sum(stage["evaluation_budget"] for stage in diagnostics["stages"])
        <= diagnostics["maximum_evaluations"]
    )


def test_native_crash_sequence_retriggers_one_stateful_body():
    events = (
        CrashEvent(0.0, 0.4, seed=1, implement=0.0, contact_spread=0.8),
        CrashEvent(0.05, 0.8, location=0.0, seed=2, implement=1.0),
    )
    sequence = render_crash_sequence(CrashFit(), 0.1, events)
    assert sequence.shape == (4800,)
    assert np.isfinite(sequence).all()
    assert np.max(np.abs(sequence[2400:])) > 0


def test_native_crash_component_taps_preserve_output():
    fit = CrashFit()
    taps = render_crash_components(fit, 0.1, strength=0.7, seed=12)
    output = render_crash(fit, 0.1, strength=0.7, seed=12)
    assert taps.shape == (4800, 5)
    assert np.array_equal(taps[:, 4], output)


def test_crash_features_retain_within_region_decay_trajectory():
    sample_rate = 48_000
    time = np.arange(5 * sample_rate) / sample_rate
    noise = np.random.default_rng(7).normal(size=time.size)
    fast = perceptual_features(AudioBuffer(noise * np.exp(-6.0 * time), sample_rate))
    slow = perceptual_features(AudioBuffer(noise * np.exp(-1.5 * time), sample_rate))
    late = [
        index
        for index, name in enumerate(fast.names)
        if name.startswith("trajectory/1.500/")
    ]
    assert 12.0 * np.mean(slow.values[late] - fast.values[late]) > 30.0


def test_crash_features_remain_finite_when_no_late_window_exists():
    samples = np.zeros(4096)
    samples[-1] = 1.0
    features = perceptual_features(AudioBuffer(samples, 48_000))
    assert np.isfinite(features.values).all()


def test_crash_features_do_not_treat_stationary_recording_noise_as_tail():
    sample_rate = 48_000
    time = np.arange(10 * sample_rate) / sample_rate
    random = np.random.default_rng(19)
    transient = random.normal(size=time.size) * np.exp(-2.5 * time)
    recording_noise = 0.01 * random.normal(size=time.size)
    features = perceptual_features(
        AudioBuffer(transient + recording_noise, sample_rate)
    )
    early = features.names.index("level/early")
    tail = features.names.index("level/tail-b")
    assert 6.0 * (features.values[early] - features.values[tail]) > 25.0


def test_crash_regrowth_gate_detects_a_delayed_energy_burst():
    sample_rate = 48_000
    time = np.arange(6 * sample_rate) / sample_rate
    natural = AudioBuffer(np.exp(-2.0 * time), sample_rate)
    burst = natural.samples.copy()
    burst[4 * sample_rate : 4 * sample_rate + sample_rate // 4] += 0.2
    assert _maximum_late_regrowth_db(natural) < 1.0
    assert _maximum_late_regrowth_db(AudioBuffer(burst, sample_rate)) > 6.0


def test_reference_residual_subtracts_an_explicit_persistent_mode():
    sample_rate = 48_000
    time = np.arange(3 * sample_rate) / sample_rate
    tau_seconds = 0.7
    samples = np.exp(-time / tau_seconds) * np.cos(2 * np.pi * 1000.0 * time)
    fit = CrashFit(
        sparse_frequency_hz=(1000.0,) + (2000.0,) * 23,
        body_decay_seconds=(np.log(1000.0) * tau_seconds,) * 8,
        sparse_decay_ratio=(1.0,) * 24,
        sparse_amplitude=(1.0,) + (0.0,) * 23,
    )
    cell = CrashFitCell("synthetic", AudioBuffer(samples, sample_rate), 1.0)
    modal, residual = reference_modal_residual(cell, fit)
    late = slice(round(0.25 * sample_rate), None)
    assert np.sqrt(np.mean(np.square(modal.samples[late]))) > 0.01
    assert np.sqrt(np.mean(np.square(residual.samples[late]))) < 1.0e-8


def test_private_cymbal_sweep_uses_crash_a_articulations():
    hits = sweep.build_sequence("private-crash-a", (32, 96), 2, 1.0, 2.0, 0.5)
    assert len(hits) == 20
    assert {(hit.articulation, hit.midi_note) for hit in hits} == {
        ("edge", 49),
        ("bow-tip", 27),
        ("bow-shank", 92),
        ("bell-tip", 93),
        ("bell-shank", 28),
    }
    midi = sweep.build_midi(hits, 0.08)
    assert midi.startswith(b"MThd")
    assert b"edge/v032/r01" in midi


def test_report_contains_only_reference_and_current_model(tmp_path):
    time = np.arange(4096) / 48000.0
    reference = AudioBuffer(np.sin(2 * np.pi * 1000 * time) * np.exp(-8 * time), 48000)
    synthesis = AudioBuffer(np.sin(2 * np.pi * 1100 * time) * np.exp(-9 * time), 48000)
    assets = tmp_path / "report-assets"
    assets.mkdir()
    (assets / "stale-reference.wav").write_bytes(b"stale")
    path = write_comparison_report(
        (ReportPair("one", reference, synthesis),),
        tmp_path / "report.html",
        "test",
        causal_fit={
            "stages": [
                {
                    "stage": "impact",
                    "end_seconds": 0.004,
                    "baseline_prefix_losses": {"0.004": 2.0},
                    "candidate_prefix_losses": {"0.004": 1.0},
                    "worst_protected_loss_ratio": 1.0,
                    "worst_current_loss_ratio": 0.5,
                    "accepted": True,
                    "absolute_quality": {
                        "passed": True,
                        "required": True,
                        "cells": [
                            {
                                "envelope_rmse_db": 1.0,
                                "spectral_rmse_db": 2.0,
                                "spectral_p95_absolute_db": 3.0,
                            }
                        ],
                    },
                }
            ],
            "completed": True,
            "blocked_stage": None,
        },
    )
    html = path.read_text(encoding="utf-8")
    assert ">Reference<" in html
    assert ">TriggerFish<" in html
    assert "old model" not in html.lower()
    assert "Reference attack, first 30 ms" in html
    assert "Reference initial decay, first 100 ms" in html
    assert "First 100 ms only" in html
    assert "First 100 ms repeated four times" in html
    assert "Cumulative causal fit" in html
    assert not (assets / "stale-reference.wav").exists()
    assert (tmp_path / "report-assets/1-one-reference.wav").is_file()
    assert (tmp_path / "report-assets/1-one-triggerfish.wav").is_file()
    assert (tmp_path / "report-assets/1-one-first-100ms-reference.wav").is_file()
    assert (tmp_path / "report-assets/1-one-first-100ms-triggerfish.wav").is_file()
    assert (tmp_path / "report-assets/1-one-first-100ms-repeat-reference.wav").is_file()
    assert (
        tmp_path / "report-assets/1-one-first-100ms-repeat-triggerfish.wav"
    ).is_file()


def test_report_spectrograms_share_one_level_reference():
    louder, quieter = _common_magnitude_db(
        np.asarray([[1.0, 0.1]]), np.asarray([[0.1, 0.01]]), -100.0
    )
    assert louder[0, 0] == 0.0
    assert quieter[0, 0] == -20.0


def test_velocity_grid_audition_gain_weights_each_velocity_once():
    quiet = AudioBuffer(np.asarray([0.01]), 48000)
    loud = AudioBuffer(np.asarray([0.1]), 48000)
    gain = velocity_grid_audition_gain(
        ((0.25, quiet), (0.25, quiet), (0.75, loud)),
        target_average_peak_dbfs=-20.0,
    )
    assert gain == pytest.approx(10.0**0.5)


def test_report_applies_one_amplifying_gain_to_both_sides(tmp_path):
    reference_samples = np.zeros(4096)
    synthesis_samples = np.zeros(4096)
    reference_samples[100] = 0.01
    synthesis_samples[100] = 0.02
    reference = AudioBuffer(reference_samples, 48000)
    synthesis = AudioBuffer(synthesis_samples, 48000)
    write_comparison_report(
        (ReportPair("quiet", reference, synthesis),),
        tmp_path / "gain-report.html",
        "gain test",
        audition_gain=10.0,
    )
    assets = tmp_path / "gain-report-assets"
    rendered_reference = read_wav(assets / "1-quiet-reference.wav", "mean")
    rendered_synthesis = read_wav(assets / "1-quiet-triggerfish.wav", "mean")
    assert np.max(np.abs(rendered_reference.samples)) == pytest.approx(0.1)
    assert np.max(np.abs(rendered_synthesis.samples)) == pytest.approx(0.2)


def test_capture_isolation_preserves_level_and_labels_cells(tmp_path):
    sample_rate = 48000
    samples = np.zeros(2 * sample_rate)
    samples[round(0.5 * sample_rate)] = 0.25
    samples[round(1.5 * sample_rate)] = 0.75
    sweep_manifest = {
        "hits": [
            {
                "index": 0,
                "onset_seconds": 0.5,
                "articulation": "edge",
                "velocity": 32,
                "repeat": 1,
            },
            {
                "index": 1,
                "onset_seconds": 1.5,
                "articulation": "bell-shank",
                "velocity": 96,
                "repeat": 2,
            },
        ]
    }
    output = tmp_path / "cells"
    cells = isolate_capture(
        AudioBuffer(samples, sample_rate), sweep_manifest, output, 0.2
    )
    manifest = output / "cells.json"
    write_cells_manifest(manifest, cells, tmp_path / "capture.wav")
    loaded = load_fit_cells(manifest)
    assert len(loaded) == 2
    assert loaded[0].location == 1.0
    assert loaded[1].location == 0.0
    assert loaded[1].hardness > loaded[0].hardness
    assert np.max(np.abs(loaded[0].reference.samples)) == np.float32(0.25)
    peak = np.argmax(np.abs(loaded[0].reference.samples))
    assert abs(peak - round(0.05 * sample_rate)) < round(0.001 * sample_rate)
    assert cells[0].cell_onset_seconds == 0.05
