"""Fixed starts and bounds shared by every perceptual-loss trial."""

from kick_control_span import coverage_layout


def experiment_starts(initial, descriptors):
    covered = coverage_layout(initial)
    covered.update(equalizer_mode=1, low_cut_hz=10, high_cut_hz=2500, colour_gain_db=0)
    starts = {
        "direct-long": dict(covered, contact_level=0.5),
        "direct-short": dict(
            covered,
            contact_level=1,
            contact_width_seconds=0.002,
            contact_noise_level=0.3,
            contact_noise_decay_seconds=0.035,
        ),
    }
    bounds = dict(
        contact_level=(0, 4),
        contact_width_seconds=(0.0005, 0.03),
        contact_noise_level=(0, 4),
        contact_noise_decay_seconds=(0.005, 0.4),
        contact_colour=(0, 1),
        thump_level=(0, 4),
        thump_pitch_hz=(20, 100),
        thump_pitch_drop_octaves=(0, 3),
        thump_pitch_fall_seconds=(0.003, 0.3),
        thump_decay_seconds=(0.08, 1),
        resonance_level=(0, 12),
        resonance_decay_seconds=(0.05, 2),
        resonance_decay_tilt=(-1, 1),
        tension_octaves=(0, 0.6),
        tension_recovery_seconds=(0.005, 0.3),
        high_cut_hz=(600, 12000),
    )
    for i in range(16):
        frequency = covered[f"resonance_frequency_{i}"]
        bounds[f"resonance_frequency_{i}"] = (
            max(20, frequency * 0.65),
            min(15000, frequency * 1.5),
        )
        if i != 1:  # one fixed prominence anchor removes a gain nullspace
            bounds[f"resonance_level_{i}"] = (-71.999, 6)
    limits = {d["key"]: d for d in descriptors}
    bounds = {
        k: (max(lo, limits[k]["minimum"]), min(hi, limits[k]["maximum"]))
        for k, (lo, hi) in bounds.items()
    }
    return starts, bounds
