"""Small-signal reference equations for the complete TB-303 filter model.

The transfer function follows Tim Stinchcombe's published analysis of the
ladder and its six surrounding AC-coupling poles/zeroes. His compact polynomial
uses an ideal 1:2 capacitor ratio; this reference substitutes the service-note
18 nF/33 nF values used by the nonlinear C++ model. The model is expected to
converge to this response at low signal levels.
"""

import numpy as np

FIRST_STAGE_SCALE = 33.0 / 18.0
CAPACITOR_SCALE = FIRST_STAGE_SCALE ** (-0.25)
STOCK_INPUT_SCALE = 0.10532968190436065
RACK_OUTPUT_SCALE = 9.494
STOCK_RESONANCE_SCALE = 0.78
HIGH_RESONANCE_MULTIPLIER = 2.0
BASS_POLE_RADIANS = 578.1
OUTPUT_KNEE_VOLTS = 8.0
OUTPUT_RAIL_VOLTS = 11.0

FORWARD_SECTIONS = (
    (578.1, 0.0),
    (97.5, 0.0),
    (20.0, 109.9),
    (4.45, 34.0),
)
OUTPUT_COUPLING_SECTIONS = ((38.5, 0.0),)
FEEDBACK_SECTIONS = (
    (578.1, 0.0),
    (97.5, 0.0),
    (38.5, 0.0),
    (20.0, 0.0),
    (7.41, 46.5),
    (4.45, 4.40),
)


def transfer(frequency, cutoff, resonance=0.0, high_resonance=False, bass=0.0):
    """Return the raw circuit's continuous-time small-signal response."""

    s = 2j * np.pi * np.asarray(frequency, dtype=float)
    omega_c = 2.0 * np.pi * cutoff
    normalized = s / omega_c
    ladder_denominator = (
        normalized**4
        + (FIRST_STAGE_SCALE + 6.0) * CAPACITOR_SCALE * normalized**3
        + (5.0 * FIRST_STAGE_SCALE + 10.0) * CAPACITOR_SCALE**2 * normalized**2
        + (6.0 * FIRST_STAGE_SCALE + 4.0) * CAPACITOR_SCALE**3 * normalized
        + 1.0
    )

    coupling_denominator = (
        (s + 97.5) * (s + 38.5) * (s + 4.45) * (s + 578.1) * (s + 20.0) * (s + 7.41)
    )
    forward_numerator = 1.06 * s**3 * (s + 109.9) * (s + 34.0) * (s + 7.41)
    feedback_numerator = 18.7 * s**4 * (s + 46.5) * (s + 4.40)
    feedback = resonance * STOCK_RESONANCE_SCALE
    if high_resonance:
        feedback *= HIGH_RESONANCE_MULTIPLIER

    response = forward_numerator / (
        ladder_denominator * coupling_denominator + feedback * feedback_numerator
    )

    bass = np.clip(bass, 0.0, 1.0)
    varied_pole = BASS_POLE_RADIANS * 10.0 ** (-bass)
    response *= (s + BASS_POLE_RADIANS) / (s + varied_pole)
    return response


def resonance_makeup(resonance, high_resonance=False):
    """Return the model's deliberately modest resonance level compensation."""

    makeup_range = 3.0 if high_resonance else 2.0
    return 1.0 + np.clip(resonance, 0.0, 1.0) * makeup_range


def analog_output_stage(voltage):
    """Return the memoryless output-compliance curve used by the C++ model."""

    values = np.asarray(voltage)
    magnitude = np.abs(values)
    headroom = OUTPUT_RAIL_VOLTS - OUTPUT_KNEE_VOLTS
    curved = OUTPUT_KNEE_VOLTS + headroom * np.tanh(
        (magnitude - OUTPUT_KNEE_VOLTS) / headroom
    )
    return np.where(magnitude <= OUTPUT_KNEE_VOLTS, values, np.copysign(curved, values))


def _cascade(input_value, states, sections, gain=1.0):
    """Evaluate a continuous analog ratio cascade and its state derivatives."""

    output = input_value
    derivatives = np.empty_like(states)
    for index, ((pole, zero), state) in enumerate(zip(sections, states)):
        derivatives[index] = pole * (output - state)
        output += (zero / pole - 1.0) * state
    return gain * output, derivatives


def render_nonlinear_reference(
    input_function,
    sample_rate,
    duration,
    cutoff,
    resonance=0.0,
    high_resonance=False,
    drive_gain=1.0,
    bass=0.0,
):
    """Render the complete continuous nonlinear model with SciPy DOP853.

    This is deliberately independent of the C++ midpoint implementation. It
    represents every analog ratio section as a continuous one-pole state and
    integrates those states together with the four nonlinear ladder states.
    """

    from scipy.integrate import solve_ivp

    forward_slice = slice(0, len(FORWARD_SECTIONS))
    feedback_slice = slice(
        forward_slice.stop, forward_slice.stop + len(FEEDBACK_SECTIONS)
    )
    ladder_slice = slice(feedback_slice.stop, feedback_slice.stop + 4)
    output_coupling_slice = slice(
        ladder_slice.stop, ladder_slice.stop + len(OUTPUT_COUPLING_SECTIONS)
    )
    bass_slice = slice(output_coupling_slice.stop, output_coupling_slice.stop + 2)
    state_size = bass_slice.stop
    feedback_amount = resonance * STOCK_RESONANCE_SCALE
    if high_resonance:
        feedback_amount *= HIGH_RESONANCE_MULTIPLIER
    varied_bass_pole = BASS_POLE_RADIANS * 10.0 ** (-np.clip(bass, 0.0, 1.0))
    bass_sections = ((varied_bass_pole, BASS_POLE_RADIANS),)
    ladder_rate = CAPACITOR_SCALE * 2.0 * np.pi * cutoff

    def derivatives(time, state):
        normalized_input = STOCK_INPUT_SCALE * drive_gain * input_function(time)
        forward, forward_derivatives = _cascade(
            normalized_input, state[forward_slice], FORWARD_SECTIONS, 1.06
        )
        ladder = state[ladder_slice]
        feedback, feedback_derivatives = _cascade(
            ladder[3], state[feedback_slice], FEEDBACK_SECTIONS, 18.7
        )
        input_junction = forward - feedback_amount * feedback
        junctions = np.tanh(
            (
                input_junction,
                ladder[0] - ladder[1],
                ladder[1] - ladder[2],
                ladder[2] - ladder[3],
                ladder[3],
            )
        )
        ladder_derivatives = ladder_rate * np.asarray(
            (
                FIRST_STAGE_SCALE * (junctions[0] - junctions[1]),
                junctions[1] - junctions[2],
                junctions[2] - junctions[3],
                junctions[3] - junctions[4],
            )
        )
        coupled_output, output_coupling_derivatives = _cascade(
            ladder[3], state[output_coupling_slice], OUTPUT_COUPLING_SECTIONS
        )
        _, bass_derivatives = _cascade(coupled_output, state[bass_slice], bass_sections)

        result = np.empty(state_size)
        result[forward_slice] = forward_derivatives
        result[feedback_slice] = feedback_derivatives
        result[ladder_slice] = ladder_derivatives
        result[output_coupling_slice] = output_coupling_derivatives
        result[bass_slice] = bass_derivatives
        return result

    sample_times = np.arange(int(sample_rate * duration)) / sample_rate
    solution = solve_ivp(
        derivatives,
        (0.0, duration),
        np.zeros(state_size),
        method="DOP853",
        t_eval=sample_times,
        rtol=2.0e-10,
        atol=2.0e-12,
        max_step=1.0 / (4.0 * sample_rate),
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    output = np.empty(sample_times.size)
    makeup = resonance_makeup(resonance, high_resonance)
    for index in range(sample_times.size):
        ladder_output = solution.y[ladder_slice, index][3]
        coupled_output, _ = _cascade(
            ladder_output,
            solution.y[output_coupling_slice, index],
            OUTPUT_COUPLING_SECTIONS,
        )
        bass_output, _ = _cascade(
            coupled_output, solution.y[bass_slice, index], bass_sections
        )
        output[index] = analog_output_stage(RACK_OUTPUT_SCALE * makeup * bass_output)
    return sample_times, output
