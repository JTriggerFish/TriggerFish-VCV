"""Independent circuit references for the ARP 4020 ADSR and board-4 AR."""

from __future__ import annotations

import numpy as np

ADSR_ATTACK_SOURCE = 1.5
SETTLING_FRACTION = 0.05
ADSR_ATTACK_CURVE = np.log(3.0)
ORDINARY_RC_CURVE = -np.log(SETTLING_FRACTION)


def normalized_exponential(phase, magnitude):
    """Endpoint-normalized RC response, with a stable linear limit."""

    phase = np.clip(np.asarray(phase, dtype=float), 0.0, 1.0)
    if abs(magnitude) < 1.0e-8:
        return phase
    return np.expm1(-magnitude * phase) / np.expm1(-magnitude)


def adsr_attack_rc(phase, source=ADSR_ATTACK_SOURCE):
    """4020 attack capacitor charging to the 10 V comparator threshold."""

    phase = np.clip(np.asarray(phase, dtype=float), 0.0, 1.0)
    magnitude = -np.log(1.0 - 1.0 / source)
    return source * (1.0 - np.exp(-magnitude * phase))


def ordinary_rc_remaining(phase):
    """Unbounded RC tail used to assess finite production endpoints."""

    phase = np.clip(np.asarray(phase, dtype=float), 0.0, 1.0)
    return np.exp(-ORDINARY_RC_CURVE * phase)


def curve_magnitude(hardware, curve, minimum, maximum):
    """Production curve mapping, independently expressed in Python."""

    curve = float(np.clip(curve, -1.0, 1.0))
    if curve < 0.0:
        return hardware * (hardware / minimum) ** curve
    return hardware * (maximum / hardware) ** curve


def production_adsr_attack(phase, curve=0.0):
    magnitude = curve_magnitude(ADSR_ATTACK_CURVE, curve, 0.1, np.log(1000.0))
    return normalized_exponential(phase, magnitude)


def production_ordinary_segment(phase, curve=0.0):
    magnitude = curve_magnitude(ORDINARY_RC_CURVE, curve, 0.25, 8.0)
    return normalized_exponential(phase, magnitude)
