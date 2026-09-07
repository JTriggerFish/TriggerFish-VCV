"""Necessary output-shape checks before a kick fit can replace the preset."""

import json
from .audio_io import read_wav
from .band_region_audit import BandRegionAudit
from .region_spectrum_audit import RegionSpectrumAudit


def heldout_seeds(saved, count=3):
    """Use distinct uint32 seeds absent from training and the primary render."""
    primary = saved["metadata"]["event"]["seed"]
    training = saved.get("training_seeds")
    if not training:
        raise ValueError("Training seeds are required for held-out validation")
    excluded = {primary, *(primary if seed is None else seed for seed in training)}
    selected = []
    seed = primary
    while len(selected) < count:
        seed = (seed + 1) & 0xFFFFFFFF
        if seed not in excluded:
            selected.append(seed)
            excluded.add(seed)
    return selected


def check_kick_candidate(directory, *, renderer=None):
    reference = read_wav(directory / "reference.wav").mono()
    candidate = read_wav(directory / "candidate.wav").mono()
    band_audit = BandRegionAudit(reference.samples, reference.sample_rate)
    shape_audit = RegionSpectrumAudit(reference.samples, reference.sample_rate)
    band = band_audit.measure(candidate.samples)
    shape = shape_audit.measure(candidate.samples)
    heldout = []
    if renderer is not None:
        saved = json.loads((directory / "search.json").read_text(encoding="utf8"))
        for seed in heldout_seeds(saved):
            audio = renderer.render(
                saved["parameters"], saved["duration_seconds"], seed
            )
            heldout.append(
                dict(
                    seed=seed,
                    band=band_audit.measure(audio),
                    shape=shape_audit.measure(audio),
                )
            )
    result = dict(
        version="kick-output-checks-v1",
        band=band,
        shape=shape,
        eligible=band["within_3db"]
        and shape["shape_guard"]
        and all(
            row["band"]["within_3db"] and row["shape"]["shape_guard"] for row in heldout
        ),
        heldout=heldout,
        heldout_checked=renderer is not None,
        listening_approved=False,
        note="Necessary engineering checks, not a perceptual acceptance score",
    )
    (directory / "quality-checks.json").write_text(
        json.dumps(result, indent=2), encoding="utf8"
    )
    return result
