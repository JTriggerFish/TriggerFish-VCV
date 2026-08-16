"""Small-signal reference equations for the ARP 4072 filter model."""

from __future__ import annotations

import numpy as np

THERMAL_VOLTAGE = 25.85e-3
STAGE_INPUT_RESISTANCE = 12_100.0
STAGE_SHUNT_RESISTANCE = 220.0
AUDIO_INPUT_RESISTANCE = 100_000.0
AUDIO_INPUT_SHUNT = 220.0
FEEDBACK_RESISTANCE = 150_000.0
FEEDBACK_INPUT_SHUNT = 2_200.0
LIMITER_TAIL_RESISTANCE = 56_000.0
POSITIVE_SUPPLY = 15.0
LIMITER_EMITTER_DROP = 0.65
OUTPUT_LEVEL_SHIFT_GAIN = 100.0 / 13.0


def circuit_values() -> dict[str, float]:
    audio_base_scale = AUDIO_INPUT_SHUNT / (AUDIO_INPUT_RESISTANCE + AUDIO_INPUT_SHUNT)
    feedback_base_scale = FEEDBACK_INPUT_SHUNT / (
        FEEDBACK_RESISTANCE + FEEDBACK_INPUT_SHUNT
    )
    tail_current = (POSITIVE_SUPPLY - LIMITER_EMITTER_DROP) / LIMITER_TAIL_RESISTANCE
    limiter_peak = tail_current * STAGE_INPUT_RESISTANCE
    limiter_slope = limiter_peak / (2.0 * THERMAL_VOLTAGE)
    return {
        "audio_base_scale": audio_base_scale,
        "feedback_base_scale": feedback_base_scale,
        "audio_base_at_5v": 5.0 * audio_base_scale,
        "feedback_base_at_5v": 5.0 * feedback_base_scale,
        "limiter_tail_current_amps": tail_current,
        "limiter_equivalent_peak_volts": limiter_peak,
        "stage_tanh_scale_per_volt": STAGE_SHUNT_RESISTANCE
        / (STAGE_INPUT_RESISTANCE * 2.0 * THERMAL_VOLTAGE),
        "output_level_shift_gain": OUTPUT_LEVEL_SHIFT_GAIN,
        "small_signal_input_gain": (
            OUTPUT_LEVEL_SHIFT_GAIN * limiter_slope * audio_base_scale
        ),
        "small_signal_feedback_gain": (
            OUTPUT_LEVEL_SHIFT_GAIN * limiter_slope * feedback_base_scale
        ),
    }


def transfer(frequency_hz, cutoff_hz, resonance=0.0):
    """Return the linearized physical-output/input transfer function."""

    frequency_hz = np.asarray(frequency_hz, dtype=float)
    s = 1j * 2.0 * np.pi * frequency_hz
    omega = 2.0 * np.pi * cutoff_hz
    one_stage = -omega / (s + omega)
    cascade = one_stage**4
    values = circuit_values()
    input_to_first_stage = (
        values["limiter_equivalent_peak_volts"]
        * values["audio_base_scale"]
        / (2.0 * THERMAL_VOLTAGE)
    )
    feedback_to_first_stage = (
        resonance
        * values["limiter_equivalent_peak_volts"]
        * values["feedback_base_scale"]
        * OUTPUT_LEVEL_SHIFT_GAIN
        / (2.0 * THERMAL_VOLTAGE)
    )
    return (
        OUTPUT_LEVEL_SHIFT_GAIN
        * cascade
        * input_to_first_stage
        / (1.0 + cascade * feedback_to_first_stage)
    )


def midpoint_transfer(frequency_hz, cutoff_hz, resonance=0.0, sample_rate=48_000.0):
    """Linearized transfer of the implementation's midpoint state update."""

    frequency_hz = np.asarray(frequency_hz, dtype=float)
    z_inverse = np.exp(-1j * 2.0 * np.pi * frequency_hz / sample_rate)
    midpoint = 0.5 * (1.0 + z_inverse)
    gamma = 2.0 * np.tan(np.pi * cutoff_hz / sample_rate)
    denominator = 1.0 - z_inverse + gamma * midpoint
    first_stage = -gamma / denominator
    later_stage = -gamma * midpoint / denominator
    cascade = first_stage * later_stage**3
    values = circuit_values()
    input_to_first_stage = (
        values["limiter_equivalent_peak_volts"]
        * values["audio_base_scale"]
        / (2.0 * THERMAL_VOLTAGE)
    )
    feedback_to_first_stage = (
        resonance
        * values["limiter_equivalent_peak_volts"]
        * values["feedback_base_scale"]
        * OUTPUT_LEVEL_SHIFT_GAIN
        / (2.0 * THERMAL_VOLTAGE)
    )
    return (
        OUTPUT_LEVEL_SHIFT_GAIN
        * cascade
        * input_to_first_stage
        / (1.0 + cascade * feedback_to_first_stage * midpoint)
    )
