"""Bounded diagnostic alternatives, not another blind local fitting pass."""

from itertools import product


def control_probes(initial):
    for key, values in {
        "contact_level": (0, 0.03, 0.1, 0.3, 1),
        "contact_width_seconds": (0.0005, 0.001, 0.003, 0.01, 0.02),
        "contact_noise_level": (0, 0.1, 0.5, 1, 3),
        "contact_noise_decay_seconds": (0.01, 0.03, 0.07, 0.15, 0.3),
        "contact_colour": (0, 0.5, 1),
        "resonance_decay_seconds": (0.15, 0.3, 0.6, 1),
        "resonance_decay_tilt": (0, 0.5, 1),
    }.items():
        for value in values:
            yield f"{key}={value}", dict(initial, **{key: value})


def coverage_layout(initial):
    """Use existing spare handles; no extra mode capacity or new DSP controls."""
    result = dict(initial)
    spare = [i for i in range(16) if initial[f"resonance_level_{i}"] <= -71.999]
    for slot, frequency in zip(
        spare, (200, 280, 390, 1100, 1550, 2200, 3100, 4300, 6000, 8000)
    ):
        result.update(
            {
                f"resonance_frequency_{slot}": frequency,
                f"resonance_level_{slot}": -24,
            }
        )
    return result


def contact_probes(initial):
    # Direct radiation, noise colour and lifetime must be varied together to
    # cross the muted-output starting point. All observations are still explicit.
    for level, decay, noise in product(
        (0.03, 0.1, 0.3), (0.01, 0.03, 0.07), (0.1, 0.5, 1)
    ):
        yield f"contact-{level}-{decay}-{noise}", dict(
            initial,
            contact_level=level,
            contact_noise_decay_seconds=decay,
            contact_noise_level=noise,
        )


def observation_probes(initial):
    # Existing ordinary low-pass observation only; no multiband/notch EQ.
    for level, decay, cutoff in product(
        (0.05, 0.15, 0.4), (0.015, 0.04, 0.08), (1200, 2200, 4000)
    ):
        yield f"lowpass-{level}-{decay}-{cutoff}", dict(
            initial,
            contact_level=level,
            contact_noise_decay_seconds=decay,
            equalizer_mode=1,
            low_cut_hz=10,
            high_cut_hz=cutoff,
            colour_gain_db=0,
        )
