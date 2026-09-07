"""Measured peaks propose handles; full renders decide whether they help."""

from .modal_fit_initialization import spectral_mode_candidates


def metallic_starts(parameters, reference, rate):
    peaks = spectral_mode_candidates(
        reference, rate, count=16, low=60, high=14000, start=0.04, end=0.5
    )
    if not peaks:
        return []
    ranked = sorted(peaks, key=lambda item: item["power_db"], reverse=True)
    starts = []
    active = [
        parameters[f"resolved_frequency_{i}"]
        for i in range(32)
        if parameters[f"resolved_level_{i}"] > -71.99
    ]
    free = [i for i in range(32) if parameters[f"resolved_level_{i}"] <= -71.99]
    missing_low = [p for p in ranked if p["frequency"] < min(active, default=400)]
    for scale in (0.0, 0.15):
        values = dict(parameters)
        for index, peak in zip(free, missing_low[:4]):
            values[f"resolved_frequency_{index}"] = peak["frequency"]
            values[f"resolved_level_{index}"] = -9.0
            values[f"resolved_turbulence_{index}"] = scale
        if free and missing_low:
            starts.append(
                (
                    f"retain layout + {min(len(free),len(missing_low),4)} lower handles, turbulence {scale}",
                    values,
                )
            )
    for count in (8, 16):
        selected = sorted(ranked[:count], key=lambda item: item["frequency"])
        # Broad high-frequency coverage is a distinct proposal, not fictitious
        # detected modes. The packet field can represent unresolved wash there.
        frequencies = sorted(
            {p["frequency"] for p in selected} | {4000.0, 6500.0, 10000.0, 14000.0}
        )
        for spread in (1.0, 4.0):
            values = dict(parameters)
            for i in range(32):
                values[f"resolved_level_{i}"] = -72.0
                values[f"resolved_turbulence_{i}"] = 1.0
                if i < len(frequencies):
                    values[f"resolved_frequency_{i}"] = frequencies[i]
                    values[f"resolved_level_{i}"] = -6.0
            values.update(
                field_packet_spread=spread,
                field_turbulence=0.65,
                field_turbulence_slope=0.3,
                field_turbulence_centre=1200.0,
                body_brightness=-6.0,
                body_excitation_centre=1200.0,
            )
            starts.append(
                (
                    f"{len(frequencies)} measured/coverage handles, spread {spread}",
                    values,
                )
            )
    return starts
