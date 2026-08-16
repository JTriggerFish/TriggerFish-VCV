import math

THERMAL_VOLTAGE = 25.85e-3
AUDIO_INPUT_RESISTANCE_OHMS = 100_000.0
AUDIO_INPUT_SHUNT_OHMS = 220.0
OUTPUT_FEEDBACK_RESISTANCE_OHMS = 56_000.0
OUTPUT_FEEDBACK_CAPACITANCE_FARADS = 100.0e-12


def circuit_values():
    audio_scale = AUDIO_INPUT_SHUNT_OHMS / (
        AUDIO_INPUT_RESISTANCE_OHMS + AUDIO_INPUT_SHUNT_OHMS
    )
    unity_current = (
        2.0 * THERMAL_VOLTAGE / (OUTPUT_FEEDBACK_RESISTANCE_OHMS * audio_scale)
    )
    return {
        "audio_input_scale": audio_scale,
        "unity_control_current_amps": unity_current,
        "output_feedback_resistance_ohms": OUTPUT_FEEDBACK_RESISTANCE_OHMS,
        "output_bandwidth_hz": 1.0
        / (
            2.0
            * math.pi
            * OUTPUT_FEEDBACK_RESISTANCE_OHMS
            * OUTPUT_FEEDBACK_CAPACITANCE_FARADS
        ),
        "linear_unity_control_volts": 10.0,
        "exponential_decibels_per_volt": 10.0,
    }


def exponential_gain(control_volts):
    return 10.0 ** ((control_volts - 10.0) / 2.0)
