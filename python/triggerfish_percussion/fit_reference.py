"""One canonical fixed-onset, fixed-gain reference extraction for all fitters."""

import numpy as np


def aligned_reference(renderer, seconds):
    if not np.isfinite(seconds) or seconds <= 0:
        raise ValueError("Reference duration must be finite and positive")
    count = round(seconds * renderer.sample_rate)
    onset_seconds = (
        renderer.metadata["reference"].get("cell", {}).get("onset_seconds", 0)
    )
    if not np.isfinite(onset_seconds) or onset_seconds < 0:
        raise ValueError("Reference onset must be finite and nonnegative")
    onset = round(onset_seconds * renderer.sample_rate)
    if onset >= len(renderer.reference):
        raise ValueError("Reference onset is outside the recording")
    samples = renderer.reference[onset : onset + count]
    return np.pad(samples, (0, count - len(samples)))
