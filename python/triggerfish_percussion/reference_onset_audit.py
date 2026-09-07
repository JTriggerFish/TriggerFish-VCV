"""Reference-only timing diagnostics; never move the audio automatically."""

import numpy as np


def audit_onset(samples, sample_rate, onset_seconds):
    samples = np.asarray(samples, dtype=float)
    onset = round(onset_seconds * sample_rate)
    if onset < 0 or onset >= len(samples) or not np.isfinite(samples).all():
        raise ValueError("Invalid reference samples or onset")
    tail = samples[onset : onset + round(0.1 * sample_rate)]
    peak = float(np.max(np.abs(tail)))
    preceding = samples[max(0, onset - round(0.02 * sample_rate)) : onset]
    noise = float(np.sqrt(np.mean(preceding**2))) if len(preceding) else None
    crossings = {}
    for db in (-50, -40, -30, -20):
        where = np.flatnonzero(np.abs(tail) > max(peak * 10 ** (db / 20), 1e-20))
        crossings[str(db)] = float(where[0] / sample_rate) if len(where) else None
    rows = []
    for a, b in ((0, 0.001), (0.001, 0.003), (0.003, 0.01), (0.01, 0.03), (0.03, 0.1)):
        window = tail[round(a * sample_rate) : round(b * sample_rate)]
        if len(window):
            rows.append(
                dict(
                    seconds=[a, b],
                    rms_dbfs=float(
                        10 * np.log10(max(float(np.mean(window**2)), 1e-20))
                    ),
                )
            )
    clear_rise = np.flatnonzero(
        np.abs(tail) > max(peak * 0.01, 6 * (noise or 0), 1e-20)
    )
    delay = float(clear_rise[0] / sample_rate) if len(clear_rise) else None
    quiet = tail[: max(0, round(((delay or 0) - 0.002) * sample_rate))]
    flat_preroll = bool(
        noise is not None
        and delay is not None
        and delay > 0.003
        and len(quiet) > 0
        and np.mean(quiet**2) < 2 * noise**2
    )
    return dict(
        declared_onset_seconds=onset_seconds,
        preceding_rms=noise,
        peak_seconds=float(np.argmax(np.abs(tail)) / sample_rate),
        relative_peak_threshold_crossings_seconds=crossings,
        rms_regions=rows,
        possible_flat_preroll=flat_preroll,
        clear_rise_offset_seconds=delay,
        automatic_alignment=False,
    )
