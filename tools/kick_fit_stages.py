"""Fit explicit kick modes and exciters. Output EQ is never a fit variable."""

STAGES = [
    (
        "thump and shared damping",
        {
            "thump_pitch_hz": (20, 100),
            "thump_pitch_drop_octaves": (0, 3),
            "thump_pitch_fall_seconds": (0.003, 0.3),
            "thump_decay_seconds": (0.08, 1),
            "thump_level": (0, 4),
            "resonance_level": (0, 12),
            "resonance_decay_seconds": (0.05, 2),
            "resonance_decay_tilt": (-1, 1),
            "tension_octaves": (0, 0.6),
            "tension_recovery_seconds": (0.005, 0.3),
        },
    ),
    (
        "contact attack and noise decay",
        {
            "contact_level": (0, 3),
            "contact_width_seconds": (0.0005, 0.03),
            "contact_colour": (0, 1),
            "contact_noise_level": (0, 4),
            "contact_noise_decay_seconds": (0.005, 0.4),
        },
    ),
]


def stages_for(parameters):
    # Fit only active centres and prominence. Spatial coupling is explicitly
    # serialized but fixed for this one-location fit; per-mode T60 does not exist.
    modes = {}
    active = [i for i in range(16) if parameters[f"resonance_level_{i}"] > -71.999]
    anchor = max(active, key=lambda i: parameters[f"resonance_level_{i}"], default=None)
    for index in range(16):
        if parameters[f"resonance_level_{index}"] <= -71.999:
            continue
        frequency = parameters[f"resonance_frequency_{index}"]
        modes[f"resonance_frequency_{index}"] = (
            max(20, frequency * 0.65),
            min(15000, frequency * 1.5),
        )
        if index != anchor:
            modes[f"resonance_level_{index}"] = (-71.999, 6)
    return [STAGES[0], ("explicit modes", modes), STAGES[1]]


def refine(search, joint_only=False, fine_only=False):
    if search.parameters["equalizer_mode"] == 2:
        raise ValueError("Disable multiband EQ before fitting the kick sources/body")
    stages = stages_for(search.parameters)
    if not joint_only:
        for name, bounds in stages:
            if bounds:
                search.stage(name, bounds, 25)
    bounds = {key: value for _, stage in stages for key, value in stage.items()}
    search.stage(
        "joint sources and modes; fixed observation",
        bounds,
        45,
        difference_step=0.001 if fine_only else 0.005,
    )
