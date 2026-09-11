"""The series fitter must not regain independent ridge controls."""

import numpy as np
import pytest

from triggerfish_percussion.structured_modal_levels import StructuredModalLevels


def parameters():
    p = {f"resolved_level_{i}": -72.0 for i in range(32)}
    for i, frequency in enumerate((125, 250, 500, 1000, 2000, 4000, 8000)):
        p[f"resolved_level_{i}"] = -20.0
        p[f"resolved_frequency_{i}"] = frequency
    return p


def test_only_two_coordinates_preserve_series_shape():
    p = parameters()
    surface = StructuredModalLevels(p)
    assert surface.matrix.shape == (7, 2)
    assert surface.levels([0, 0]) == pytest.approx(surface.initial)
    assert surface.levels([2, 1]) == pytest.approx(np.arange(-21, -14))
    assert np.diff(surface.levels([2, 1]), n=2) == pytest.approx(np.zeros(5))
    assert p == parameters()  # construction/evaluation never mutate a snapshot


def test_an_isolated_ridge_cannot_be_independently_adjusted():
    surface = StructuredModalLevels(parameters())
    spike = np.zeros(7)
    spike[3] = 12
    controls = np.linalg.lstsq(surface.matrix, spike, rcond=None)[0]
    assert np.linalg.norm(surface.matrix @ controls - spike) > 10


def test_disabled_modes_stay_out_and_controls_do_not_clip():
    surface = StructuredModalLevels(parameters())
    assert len(surface.keys) == 7
    assert np.max(surface.levels([40, 0])) == 20
    with pytest.raises(ValueError):
        surface.levels(np.zeros(7))
    with pytest.raises(ValueError):
        surface.levels([0, np.nan])


def test_empty_series_rejected():
    with pytest.raises(ValueError):
        StructuredModalLevels({f"resolved_level_{i}": -72 for i in range(32)})
