"""Shared-coordinate search stages; no hidden runtime parameters or publication."""

import numpy as np
from scipy.optimize import minimize


def search_metadata(args):
    """Name the actual fitting coordinates in the saved audit trail."""
    if args.low_decay:
        name = "low endpoint refinement"
        keys = ["low T60 multiplier"]
    elif args.locked_texture:
        name = "texture-locked tuning/decay"
        keys = ["low T60", "high T60", "pitch scale", "upper stretch"]
    elif args.upper_balance:
        name = "upper prominence"
        keys = ["upper prominence dB", "upper slope endpoint Hz"]
    elif args.sparse_decay:
        name = "fine shared decay" if args.fine_decay else "sparse shared decay"
        keys = ["knot frequency", "knot T60 multiplier", "upper T60 multiplier"]
        if not args.fine_decay:
            keys = ["energy sensitivity", "bloom rate multiplier"] + keys
    else:
        name = "shared series"
        keys = [
            "pitch scale",
            "upper stretch",
            "low T60",
            "high T60",
            "bloom rate",
            "excitation tilt",
            "excitation centre",
            "bass shelf",
            "mid shelf",
        ]
    return dict(search_stage=name, search_coordinates=keys)


def shared_edit(base, values):
    """Global pitch/stretch and broad prominence shelves; no individual ridge search."""
    values = np.asarray(values, dtype=float)
    if values.shape != (9,) or not np.isfinite(values).all():
        raise ValueError("Expected nine finite shared coordinates")
    p = dict(base)
    scale, stretch, low, high, bloom, tilt, centre, bass, middle = values
    p.update(
        body_decay_seconds_0=low,
        body_decay_seconds_7=high,
        bloom_rate=bloom,
        body_brightness=tilt,
        body_excitation_centre=centre,
    )
    for i in range(32):
        if base[f"resolved_level_{i}"] <= -71.99:
            continue
        f = base[f"resolved_frequency_{i}"]
        # Concentrate stretching above 500 Hz, with zero offset at 120 Hz.
        bend = np.log2(1 + (f / 500) ** 2) - np.log2(1 + (120 / 500) ** 2)
        p[f"resolved_frequency_{i}"] = float(
            np.clip(f * scale * 2 ** (stretch * bend), 1, 15000)
        )
        db = np.interp(
            np.log(f),
            np.log([100, 300, 1000, 3000, 15000]),
            [bass, bass, middle, middle, 0],
        )
        if db != 0:
            p[f"resolved_level_{i}"] = float(
                np.clip(base[f"resolved_level_{i}"] + db, -72, 6)
            )
    return p


def texture_cases(base):
    yield "unchanged", base
    for depth in (0.3, 0.7, 1.2):
        for rate in (20, 60, 160):
            yield f"shimmer {depth}/{rate}", dict(
                base,
                output_eq_enabled=0,
                field_phase_bandwidth=0,
                field_motion_depth=depth,
                field_motion_rate=rate,
                field_motion_sharing=0.15,
                field_satellite_density=1,
                field_beat_depth=0.15,
            )
    yield "EQ bypass only", dict(base, output_eq_enabled=0)
    yield "quiet paired rings", dict(
        base,
        output_eq_enabled=0,
        field_phase_bandwidth=0,
        field_beat_depth=0.1,
        field_satellite_density=1,
    )


def damping_cases(base, fine=False):
    """One shared knot only, tested after geometry and two-endpoint fitting."""
    if any(v >= 0.5 for k, v in base.items() if k.startswith("body_decay_active_")):
        raise ValueError("Sparse damping trial expects a two-endpoint starting curve")
    erb = lambda f: 21.4 * np.log10(1 + 0.00437 * f)
    seen = set()
    for frequency in (150, 600, 3000):
        amount = (erb(frequency) - erb(40)) / (erb(15000) - erb(40))
        current = np.exp(
            (1 - amount) * np.log(base["body_decay_seconds_0"])
            + amount * np.log(base["body_decay_seconds_7"])
        )
        for scale in ((0.9, 1, 1.1) if fine else (0.6, 1.4, 2)):
            for upper in ((0.7, 1) if fine else (0.5, 1)):
                seconds = float(np.clip(current * scale, 0.02, 30))
                key = (frequency, seconds, upper)
                if key in seen:
                    continue
                seen.add(key)
                yield str(key), dict(
                    base,
                    body_decay_active_1=1,
                    body_decay_frequency_1=frequency,
                    body_decay_seconds_1=seconds,
                    body_decay_seconds_7=base["body_decay_seconds_7"] * upper,
                )


def fit_sparse_decay(initial, rows, evaluate, fine=False):
    """Try existing transfer controls before spending one extra damping knot."""
    for power in (() if fine else (0.7, 1, 1.4)):
        for scale in (0.6, 1):
            evaluate(
                dict(
                    initial,
                    bloom_energy_sensitivity=power,
                    bloom_rate=initial["bloom_rate"] * scale,
                ),
                f"late transfer {power}/{scale}",
            )
    base = min(rows, key=lambda r: r["score"])["parameters"]
    for name, parameters in damping_cases(base, fine):
        evaluate(parameters, "one shared decay knot " + name)


def fit_shared_series(initial, rows, evaluate, budget):
    """Screen texture, then search nine scaled coordinates without ridge chasing."""
    for name, parameters in texture_cases(initial):
        if name != "unchanged":
            evaluate(parameters, name)
    base = min(
        (r for r in rows if r["parameters"]["output_eq_enabled"] == 0),
        key=lambda r: r["score"],
    )["parameters"]
    bounds = [
        (0.8, 1.15),
        (-0.075, 0.075),
        (5, 30),
        (0.15, 3),
        (2, 12),
        (-18, 0),
        (800, 4000),
        (-6, 6),
        (-6, 6),
    ]
    start = [
        1,
        0,
        base["body_decay_seconds_0"],
        base["body_decay_seconds_7"],
        base["bloom_rate"],
        base["body_brightness"],
        base["body_excitation_centre"],
        0,
        0,
    ]
    lower, upper = np.asarray(bounds).T

    def score(x):
        return evaluate(
            shared_edit(base, lower + np.asarray(x) * (upper - lower)),
            "shared series/dynamics",
        )

    minimize(
        score,
        np.clip((np.asarray(start) - lower) / (upper - lower), 0, 1),
        method="Powell",
        bounds=[(0, 1)] * len(bounds),
        options=dict(maxfev=budget, xtol=0.015, ftol=0.002),
    )


def fit_upper_balance(initial, evaluate):
    """Broad painted prominence slope, not output EQ or individual ridge fitting."""
    for end in (4000, 6000):
        for level in (-2, -4, -6):
            p = dict(initial)
            for index in range(32):
                key = f"resolved_level_{index}"
                if p[key] <= -71.99:
                    continue
                frequency = p[f"resolved_frequency_{index}"]
                weight = np.clip(np.log(frequency / 1500) / np.log(end / 1500), 0, 1)
                p[key] = float(max(-72, p[key] + level * weight))
            evaluate(p, f"upper prominence {level} dB by {end} Hz")


def tuning_decay_edit(base, values):
    """Four shared coordinates; retain all user texture, gain and transport settings."""
    low, high, pitch, stretch = values
    return shared_edit(
        base,
        [
            pitch,
            stretch,
            low,
            high,
            base["bloom_rate"],
            base["body_brightness"],
            base["body_excitation_centre"],
            0,
            0,
        ],
    )


def fit_locked_texture(initial, evaluate, budget):
    """Refine damping first, then tuning, without undoing a listening-approved texture."""
    lower = np.array([12, 0.25, 0.96, -0.03])
    upper = np.array([30, 3, 1.04, 0.03])
    start = np.array(
        [initial["body_decay_seconds_0"], initial["body_decay_seconds_7"], 1, 0]
    )

    def score(x):
        return evaluate(
            tuning_decay_edit(initial, lower + x * (upper - lower)),
            "texture-locked tuning/decay",
        )

    minimize(
        score,
        np.clip((start - lower) / (upper - lower), 0, 1),
        method="Powell",
        bounds=[(0, 1)] * 4,
        options=dict(maxfev=budget, xtol=0.01, ftol=0.002),
    )
