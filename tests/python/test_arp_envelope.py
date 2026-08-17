import numpy as np

import _triggerfish_dsp as dsp
from reference_arp_envelope import (
    adsr_attack_rc,
    ordinary_rc_remaining,
    production_adsr_attack,
    production_ordinary_segment,
)

SAMPLE_RATE = 48_000.0


def _render_gate(
    duration,
    *,
    attack,
    decay=0.2,
    sustain=0.4,
    release=0.3,
    curve=0.0,
    mode=0,
    trigger=None,
    auto_gate_trigger=True,
):
    count = int(round(duration * SAMPLE_RATE))
    gate = np.full(count, 10.0)
    if trigger is None:
        trigger = np.zeros(count)
    return dsp.arp_envelope(
        gate,
        trigger,
        attack,
        decay,
        sustain,
        release,
        curve,
        mode,
        auto_gate_trigger,
        SAMPLE_RATE,
    )


def test_default_adsr_attack_matches_original_15v_threshold_reference():
    attack = 0.123456
    output = _render_gate(attack, attack=attack)
    phase = np.arange(1, len(output) + 1) / (attack * SAMPLE_RATE)
    expected = adsr_attack_rc(phase)
    assert np.max(np.abs(output - expected)) < 2.0e-12
    assert output[-1] == 1.0


def test_default_falling_curve_tracks_the_rc_to_its_five_percent_endpoint():
    phase = np.linspace(0.0, 1.0, 1001)
    production_remaining = 1.0 - production_ordinary_segment(phase)
    circuit_remaining = ordinary_rc_remaining(phase)
    error = production_remaining - circuit_remaining
    assert np.sqrt(np.mean(error**2)) < 0.039
    assert np.max(np.abs(error)) <= 0.05 + 1.0e-12


def test_curve_changes_shape_without_changing_attack_duration():
    attacks = [_render_gate(0.1, attack=0.1, curve=curve) for curve in (-1.0, 0.0, 1.0)]
    quarter = int(0.025 * SAMPLE_RATE) - 1
    assert attacks[0][quarter] < attacks[1][quarter] < attacks[2][quarter]
    assert all(output[-1] == 1.0 for output in attacks)
    for curve, output in zip((-1.0, 0.0, 1.0), attacks, strict=True):
        phase = np.arange(1, len(output) + 1) / (0.1 * SAMPLE_RATE)
        assert np.max(np.abs(output - production_adsr_attack(phase, curve))) < 2e-12


def test_ar_uses_ordinary_rc_attack_and_ignores_trigger():
    trigger = np.zeros(int(0.1 * SAMPLE_RATE))
    trigger[1200:1202] = 10.0
    output = _render_gate(0.1, attack=0.1, mode=2, trigger=trigger)
    phase = np.arange(1, len(output) + 1) / (0.1 * SAMPLE_RATE)
    expected = production_ordinary_segment(phase)
    assert np.max(np.abs(output - expected)) < 2.0e-12


def test_patched_trigger_enables_strict_gate_and_trigger_behavior():
    count = int(0.3 * SAMPLE_RATE)
    gate = np.full(count, 10.0)
    trigger = np.zeros(count)
    trigger_sample = int(0.1 * SAMPLE_RATE)
    trigger[trigger_sample : trigger_sample + 2] = 10.0
    output = dsp.arp_envelope(
        gate,
        trigger,
        0.05,
        0.1,
        0.4,
        0.1,
        0.0,
        0,
        False,
        SAMPLE_RATE,
    )
    assert 0.39 < output[trigger_sample - 1] <= 0.4
    assert output[trigger_sample] > output[trigger_sample - 1]
    assert output[-1] == 0.4


def test_ad_runs_attack_then_decay_without_waiting_for_gate_fall():
    attack = 0.05
    decay = 0.1
    output = _render_gate(
        0.3,
        attack=attack,
        decay=decay,
        sustain=0.8,
        release=0.25,
        mode=1,
    )
    peak = int(round(attack * SAMPLE_RATE)) - 1
    end = int(round((attack + decay) * SAMPLE_RATE)) - 1
    assert output[peak] == 1.0
    assert output[peak + 1] < 1.0
    assert output[end] == 0.0
    assert np.all(output[end:] == 0.0)


def test_five_millisecond_attack_remains_available():
    output = _render_gate(0.005, attack=0.005)
    assert output[-1] == 1.0
