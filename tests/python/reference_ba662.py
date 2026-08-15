"""Published reduced equations and calibration bounds for the BA662 study."""

import numpy as np

THERMAL_VOLTAGE = 0.02585
CURRENT_TRANSFER_EFFICIENCY = 0.85
RACK_TO_DIFFERENTIAL_SCALE = np.sqrt(2.0) * 1.0e-3
BASE_CONTROL_CURRENT = 20.0e-6
ACCENT_CONTROL_CURRENT = 20.0e-6
PHYSICAL_OUTPUT_LOAD_OHMS = 220_000.0
RACK_OUTPUT_CALIBRATION = 9.818181818181818
RACK_TRANSIMPEDANCE_OHMS = PHYSICAL_OUTPUT_LOAD_OHMS * RACK_OUTPUT_CALIBRATION


def c38_transfer(
    frequency_hz,
    *,
    capacitance_farads=1.0e-6,
    load_ohms=50_000.0,
    source_ohms=0.0,
    esr_ohms=0.0,
    leakage_ohms=np.inf,
):
    """Analog C38/load transfer with optional capacitor sensitivity terms."""

    frequency_hz = np.asarray(frequency_hz, dtype=float)
    leakage_admittance = 0.0 if np.isinf(leakage_ohms) else 1.0 / leakage_ohms
    capacitor_impedance = esr_ohms + 1.0 / (
        2j * np.pi * frequency_hz * capacitance_farads + leakage_admittance
    )
    return load_ohms / (source_ohms + capacitor_impedance + load_ohms)


def differential_pair_current(
    differential_voltage,
    control_current,
    *,
    efficiency=CURRENT_TRANSFER_EFFICIENCY,
    thermal_voltage=THERMAL_VOLTAGE,
):
    """Exact matched-BJT-pair output current before mirror imperfections."""

    differential_voltage = np.asarray(differential_voltage, dtype=float)
    control_current = np.clip(np.asarray(control_current, dtype=float), 0.0, 2.0e-3)
    return (
        efficiency
        * control_current
        * np.tanh(differential_voltage / (2.0 * thermal_voltage))
    )


def small_signal_gm(
    control_current,
    *,
    efficiency=CURRENT_TRANSFER_EFFICIENCY,
    thermal_voltage=THERMAL_VOLTAGE,
):
    return efficiency * np.asarray(control_current) / (2.0 * thermal_voltage)


def static_tb303_vca(audio_rack_volts, base_control, accent_control):
    """Static wrapper transfer, excluding only the sub-audio C38 high-pass."""

    control_current = BASE_CONTROL_CURRENT * np.clip(
        base_control, 0.0, 1.0
    ) + ACCENT_CONTROL_CURRENT * np.clip(accent_control, 0.0, 1.0)
    current = differential_pair_current(
        RACK_TO_DIFFERENTIAL_SCALE * np.asarray(audio_rack_volts),
        control_current,
    )
    return RACK_TRANSIMPEDANCE_OHMS * current
