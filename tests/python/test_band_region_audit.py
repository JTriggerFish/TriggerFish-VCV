"""Independent region checks must expose spectral holes and excess late energy."""

import numpy as np
from triggerfish_percussion.band_region_audit import BandRegionAudit

RATE = 24000
TIME = np.arange(round(1.2 * RATE)) / RATE


def tone(frequency, decay=0.15):
    return 0.1 * np.sin(2 * np.pi * frequency * TIME) * 10 ** (-3 * TIME / decay)


def test_identity_and_missing_upper_attack():
    reference = tone(60, 0.4) + tone(1700)
    audit = BandRegionAudit(reference, RATE)
    assert audit.measure(reference)["within_3db"]
    assert not audit.measure(tone(60, 0.4))["within_3db"]


def test_sub_khz_resonance_cannot_replace_upper_attack():
    reference = tone(1700)
    result = BandRegionAudit(reference, RATE).measure(tone(800))
    assert result["error_db"][4][0] < -12


def test_silent_reference_rejects_audible_candidate():
    audit = BandRegionAudit(np.zeros_like(TIME), RATE)
    assert audit.measure(np.zeros_like(TIME))["within_3db"]
    assert not audit.measure(tone(100))["within_3db"]


def test_late_energy_is_not_excluded_by_reference_silence():
    reference = tone(100)
    late = np.where(TIME > 0.7, 0.05 * np.sin(2 * np.pi * 100 * TIME), 0)
    result = BandRegionAudit(reference, RATE).measure(reference + late)
    assert not result["within_3db"]
    assert result["error_db"][0][4] > 20
