"""Small identifiable parameter blocks; events, velocity response and routing fixed."""


def snare_stages(parameters):
    return [
        (
            "body pitch and loss",
            dict(
                fundamental_hz=(100, 320),
                decay_seconds=(0.05, 0.8),
                decay_tilt=(-1, 1),
                inharmonicity=(0, 1),
                body_brightness=(0, 1),
                tension_octaves=(0, 0.3),
            ),
            range(5),
        ),
        (
            "wire spectrum and decay",
            dict(
                wire_level=(0.05, 5),
                wire_decay_seconds=(0.015, 0.8),
                wire_decay_tilt=(-0.5, 1),
                wire_release_seconds=(0.005, 0.19),
                wire_minimum_hz=(100, 2500),
                wire_maximum_hz=(4000, 22000),
                wire_brightness=(0, 1),
                ring_level=(0, 1),
            ),
            range(5),
        ),
        (
            "contact and body balance",
            dict(
                contact_direct_level=(0, 1.5),
                contact_body_level=(0.05, 2),
                contact_duration_seconds=(0.0003, 0.012),
                contact_brightness=(0, 1),
                contact_noise_level=(0, 2),
                contact_noise_decay_seconds=(0.001, 0.06),
            ),
            range(5),
        ),
        (
            "ring and output bandwidth",
            dict(
                ring_frequency_hz=(250, 1200),
                ring_decay_seconds=(0.05, 1.5),
                ring_level=(0, 0.7),
                high_cut_hz=(6000, 22000),
            ),
            range(5),
        ),
        (
            "wire noise balance (modal mix anchored)",
            dict(
                wire_noise_mix=(0.05, 2),
                wire_release_seconds=(0.005, 0.19),
                wire_motion_highpass_hz=(30, 1000),
                decay_seconds=(0.05, 0.8),
                decay_tilt=(-1, 1),
            ),
            range(5),
        ),
    ]


def metallic_stages(parameters):
    modes = [i for i in range(32) if parameters[f"resolved_level_{i}"] > -71.99]
    # Local frequency polish targets salient ringing, not every dense-wash
    # handle. Moving broad packets by 3% adds weak, seed-sensitive directions.
    frequency_modes = sorted(
        modes, key=lambda i: parameters[f"resolved_level_{i}"], reverse=True
    )
    frequency_modes = [
        i for i in frequency_modes if parameters[f"resolved_frequency_{i}"] < 3000
    ][:8]
    return [
        (
            "two endpoint damping and energy travel",
            dict(
                body_decay_seconds_0=(0.1, 30),
                body_decay_seconds_7=(0.02, 15),
                bloom_rate=(0, 16),
                body_brightness=(-72, 24),
                body_excitation_centre=(80, 10000),
            ),
            range(5),
        ),
        (
            "modal observation spectrum",
            {f"resolved_level_{i}": (-45, 6) for i in modes},
            range(5),
        ),
        (
            "density and correlation",
            dict(
                field_turbulence=(0, 1),
                field_turbulence_slope=(-1, 1),
                field_packet_spread=(0.05, 10),
                field_phase_bandwidth=(0, 4),
            ),
            range(5),
        ),
        (
            "local mode frequencies",
            {
                f"resolved_frequency_{i}": (
                    max(40, parameters[f"resolved_frequency_{i}"] * 0.97),
                    min(15000, parameters[f"resolved_frequency_{i}"] * 1.03),
                )
                for i in frequency_modes
            },
            range(5),
        ),
        (
            "contact and radiation",
            dict(
                direct_gain=(0, 1),
                impact_width=(0.25, 4),
                impact_tone_noise=(0, 1),
                impact_noise_tilt=(-18, 18),
                output_high_cut=(3000, 22000),
            ),
            range(5),
        ),
    ]
