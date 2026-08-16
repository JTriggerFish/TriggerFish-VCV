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
OUTPUT_KNEE_VOLTS = 10.0
OUTPUT_RAIL_VOLTS = 13.5


def stage_base_resistance() -> float:
    """The exact shunt resistance after loading by both 12.1 kohm paths."""

    return 1.0 / (1.0 / STAGE_SHUNT_RESISTANCE + 2.0 / STAGE_INPUT_RESISTANCE)


def limiter_differential_resistance() -> float:
    """Mean collector load seen by the limiter's differential current."""

    return 0.5 * (STAGE_SHUNT_RESISTANCE + stage_base_resistance())


def circuit_values() -> dict[str, float]:
    audio_base_scale = AUDIO_INPUT_SHUNT / (AUDIO_INPUT_RESISTANCE + AUDIO_INPUT_SHUNT)
    feedback_base_scale = FEEDBACK_INPUT_SHUNT / (
        FEEDBACK_RESISTANCE + FEEDBACK_INPUT_SHUNT
    )
    tail_current = (POSITIVE_SUPPLY - LIMITER_EMITTER_DROP) / LIMITER_TAIL_RESISTANCE
    nominal_limiter_peak = (
        tail_current
        * STAGE_INPUT_RESISTANCE
        * limiter_differential_resistance()
        / stage_base_resistance()
    )
    unity_limiter_peak = (
        2.0 * THERMAL_VOLTAGE / (OUTPUT_LEVEL_SHIFT_GAIN * audio_base_scale)
    )
    limiter_calibration = unity_limiter_peak / nominal_limiter_peak
    limiter_peak = nominal_limiter_peak * limiter_calibration
    limiter_slope = limiter_peak / (2.0 * THERMAL_VOLTAGE)
    return {
        "audio_base_scale": audio_base_scale,
        "feedback_base_scale": feedback_base_scale,
        "audio_base_at_5v": 5.0 * audio_base_scale,
        "feedback_base_at_5v": 5.0 * feedback_base_scale,
        "limiter_tail_current_amps": tail_current,
        "limiter_equivalent_peak_volts": limiter_peak,
        "nominal_limiter_equivalent_peak_volts": nominal_limiter_peak,
        "limiter_gain_calibration": limiter_calibration,
        "limiter_differential_resistance_ohms": limiter_differential_resistance(),
        "stage_saturation_coefficient_per_volt": stage_base_resistance()
        / (STAGE_INPUT_RESISTANCE * 2.0 * THERMAL_VOLTAGE),
        "stage_base_resistance_ohms": stage_base_resistance(),
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


def output_compliance(voltage):
    """Final level-shifter compliance used by the production filter."""

    voltage = np.asarray(voltage, dtype=float)
    magnitude = np.abs(voltage)
    headroom = OUTPUT_RAIL_VOLTS - OUTPUT_KNEE_VOLTS
    limited = np.sign(voltage) * (
        OUTPUT_KNEE_VOLTS
        + headroom * np.tanh((magnitude - OUTPUT_KNEE_VOLTS) / headroom)
    )
    return np.where(magnitude <= OUTPUT_KNEE_VOLTS, voltage, limited)


def continuous_rhs(
    _time,
    state,
    input_value,
    cutoff_hz,
    resonance,
    drive_gain=1.0,
    exact_loading=True,
):
    """Continuous four-state 4072 equations used by the DOP853 reference."""

    resistance = stage_base_resistance() if exact_loading else STAGE_SHUNT_RESISTANCE
    stage_saturation_coefficient = resistance / (
        STAGE_INPUT_RESISTANCE * 2.0 * THERMAL_VOLTAGE
    )
    values = circuit_values()
    limiter_voltage = (
        drive_gain * input_value * values["audio_base_scale"]
        - resonance * values["feedback_base_scale"] * OUTPUT_LEVEL_SHIFT_GAIN * state[3]
    )
    limiter_output = values["limiter_equivalent_peak_volts"] * np.tanh(
        limiter_voltage / (2.0 * THERMAL_VOLTAGE)
    )
    stage_inputs = np.asarray(
        (
            limiter_output + state[0],
            state[0] + state[1],
            state[1] + state[2],
            state[2] + state[3],
        )
    )
    return -(2.0 * np.pi * cutoff_hz / stage_saturation_coefficient) * np.tanh(
        stage_saturation_coefficient * stage_inputs
    )


def render_continuous_sine(
    frequency_hz,
    amplitude,
    cutoff_hz,
    resonance,
    *,
    drive_gain=1.0,
    duration=0.25,
    sample_rate=48_000.0,
    exact_loading=True,
    rtol=2.0e-11,
    atol=2.0e-12,
    max_step_divisor=8.0,
):
    """Render the continuous nonlinear circuit reduction with SciPy DOP853."""

    from scipy.integrate import solve_ivp

    sample_times = np.arange(int(sample_rate * duration)) / sample_rate

    def derivatives(time, state):
        input_value = amplitude * np.sin(2.0 * np.pi * frequency_hz * time)
        return continuous_rhs(
            time,
            state,
            input_value,
            cutoff_hz,
            resonance,
            drive_gain,
            exact_loading,
        )

    solution = solve_ivp(
        derivatives,
        (0.0, duration),
        np.zeros(4),
        method="DOP853",
        t_eval=sample_times,
        rtol=rtol,
        atol=atol,
        max_step=1.0 / (max_step_divisor * sample_rate),
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return sample_times, output_compliance(OUTPUT_LEVEL_SHIFT_GAIN * solution.y[3])


def render_continuous_self_oscillation(
    cutoff_hz,
    resonance,
    *,
    duration=2.0,
    sample_rate=48_000.0,
    exact_loading=True,
    rtol=2.0e-11,
    atol=2.0e-12,
    max_step_divisor=8.0,
):
    """Excite the autonomous continuous model with a microscopic state impulse."""

    from scipy.integrate import solve_ivp

    sample_times = np.arange(int(sample_rate * duration)) / sample_rate

    def derivatives(time, state):
        return continuous_rhs(
            time,
            state,
            0.0,
            cutoff_hz,
            resonance,
            1.0,
            exact_loading,
        )

    initial = np.asarray((1.0e-9, 0.0, 0.0, 0.0))
    solution = solve_ivp(
        derivatives,
        (0.0, duration),
        initial,
        method="DOP853",
        t_eval=sample_times,
        rtol=rtol,
        atol=atol,
        max_step=1.0 / (max_step_divisor * sample_rate),
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return sample_times, output_compliance(OUTPUT_LEVEL_SHIFT_GAIN * solution.y[3])


def harmonic_amplitudes(signal, frequency_hz, sample_rate, count=9, cycles=100):
    """Project the final whole cycles onto integer harmonics."""

    samples_per_cycle = sample_rate / frequency_hz
    selected_count = min(len(signal), int(round(cycles * samples_per_cycle)))
    selected = np.asarray(signal[-selected_count:], dtype=float)
    time = np.arange(selected_count) / sample_rate
    return np.asarray(
        [
            2.0
            * abs(
                np.mean(selected * np.exp(-2j * np.pi * harmonic * frequency_hz * time))
            )
            for harmonic in range(1, count + 1)
        ]
    )


def positive_zero_crossing_frequency(signal, sample_rate):
    """Estimate a steady oscillator frequency with interpolated rising crossings."""

    signal = np.asarray(signal, dtype=float)
    indices = np.flatnonzero((signal[:-1] <= 0.0) & (signal[1:] > 0.0))
    if indices.size < 3:
        raise ValueError("at least three rising zero crossings are required")
    fractions = -signal[indices] / (signal[indices + 1] - signal[indices])
    crossings = (indices + fractions) / sample_rate
    return 1.0 / np.mean(np.diff(crossings))
