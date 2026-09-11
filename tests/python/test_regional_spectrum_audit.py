"""Known missing-tone counterexamples, not just identity-distance checks."""

import numpy as np
import pytest

from triggerfish_percussion.layered_band_loss import LayeredBandLoss
from triggerfish_percussion.regional_spectrum_audit import RegionalSpectrumAudit


@pytest.mark.parametrize("frequencies", [(120, 220), (300, 440), (1000, 1150)])
def test_broadband_improvement_can_remove_a_reference_tone(frequencies):
    rate = 32000
    t = np.arange(6 * rate) / rate
    envelope = np.exp(-t / 0.6)
    low = 0.05 * envelope * np.sin(2 * np.pi * frequencies[0] * t)
    high = 0.1 * envelope * np.sin(2 * np.pi * frequencies[1] * t)
    reference = low + high
    incoming = 3 * low + high
    # Same total sinusoidal power as reference, but the lower tone is gone.
    missing = np.sqrt(1.25) * high
    bands = LayeredBandLoss(reference, rate, audibility=True)
    audit = RegionalSpectrumAudit(reference, rate, [(0.05, 0.2), (0.2, 0.5)])
    assert bands.score(missing) < bands.score(incoming)
    assert audit.measure(incoming)["max_deficit_db"] < 1
    assert audit.measure(missing)["max_deficit_db"] > 20
    assert audit.measure(reference)["max_deficit_db"] == 0


def test_fixed_reference_mask_and_floor_cannot_follow_candidate_gain():
    rate = 32000
    t = np.arange(rate) / rate
    reference = np.sin(2 * np.pi * 170 * t) * np.exp(-t)
    audit = RegionalSpectrumAudit(reference, rate, [(0.1, 0.5)])
    before = audit.active.copy(), audit.floor.copy()
    assert audit.measure(reference / 10)["max_deficit_db"] == pytest.approx(20)
    assert audit.measure(reference * 10)["max_excess_db"] == pytest.approx(20)
    assert np.array_equal(audit.active, before[0])
    assert np.array_equal(audit.floor, before[1])
    with pytest.raises(ValueError, match="finite"):
        audit.measure(reference * np.nan)
    with pytest.raises(ValueError, match="finite"):
        audit.measure(reference[:-1])


def test_silence_is_finite_and_invalid_regions_are_rejected():
    reference = np.zeros(32000)
    audit = RegionalSpectrumAudit(reference, 32000, [(0, 0.1)])
    assert audit.measure(reference)["max_deficit_db"] == 0
    for regions in ([], [(0.1, 0)], [(0, 2)], [(0, 0.001)], [(0, np.nan)]):
        with pytest.raises(ValueError):
            RegionalSpectrumAudit(reference, 32000, regions)
    with pytest.raises(ValueError):
        RegionalSpectrumAudit(0.0, 32000, [(0, 0.1)])
    RegionalSpectrumAudit(reference, 32000, [(0.25, 0.3)])
