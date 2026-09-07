"""Spectral shape checks must reject shifted modes even at matching gain."""

import numpy as np
from triggerfish_percussion.region_spectrum_audit import RegionSpectrumAudit
from triggerfish_percussion.region_fit_loss import RegionFitLoss
from triggerfish_percussion.band_region_audit import BandRegionAudit


def test_shifted_mode_is_not_an_identity():
    rate = 48000
    t = np.arange(round(1.2 * rate)) / rate
    reference = 0.1 * np.sin(2 * np.pi * 1700 * t) * np.exp(-40 * t)
    wrong = 0.1 * np.sin(2 * np.pi * 800 * t) * np.exp(-40 * t)
    audit = RegionSpectrumAudit(reference, rate)
    assert audit.measure(reference)["shape_guard"]
    assert not audit.measure(wrong)["shape_guard"]


def test_common_sample_rates_and_silence():
    for rate in (24000, 44100, 48000, 96000):
        silence = np.zeros(round(1.2 * rate))
        audit = RegionSpectrumAudit(silence, rate)
        assert np.isfinite(audit.target).all()
        assert audit.measure(silence)["shape_guard"]
        assert not audit.measure(np.ones_like(silence) * 0.1)["shape_guard"]


def test_region_fit_keeps_absolute_gain():
    rate = 24000
    t = np.arange(round(1.2 * rate)) / rate
    reference = (
        np.sin(2 * np.pi * 80 * t) + 0.1 * np.sin(2 * np.pi * 1700 * t)
    ) * np.exp(-30 * t)
    loss = RegionFitLoss(BandRegionAudit(reference, rate), reference, rate)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(reference * 0.5)) > 1


def test_region_fit_penalizes_new_energy_in_quiet_cells():
    rate = 24000
    t = np.arange(round(1.2 * rate)) / rate
    reference = 0.1 * np.sin(2 * np.pi * 80 * t) * np.exp(-30 * t)
    loss = RegionFitLoss(BandRegionAudit(reference, rate), reference, rate)
    wrong = reference + np.where(t > 0.7, 0.03 * np.sin(2 * np.pi * 4000 * t), 0)
    assert np.linalg.norm(loss.residual(reference)) == 0
    assert np.linalg.norm(loss.residual(wrong)) > 10
