from types import SimpleNamespace
import numpy as np
import pytest

from triggerfish_percussion.fit_reference import aligned_reference


def test_fixed_onset_and_padding_without_gain_adjustment():
    renderer = SimpleNamespace(
        sample_rate=10,
        reference=np.arange(12, dtype=float),
        metadata={"reference": {"cell": {"onset_seconds": 0.2}}},
    )
    np.testing.assert_array_equal(aligned_reference(renderer, 0.5), [2, 3, 4, 5, 6])
    np.testing.assert_array_equal(
        aligned_reference(renderer, 1.2), [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0]
    )
    renderer.metadata["reference"]["cell"]["onset_seconds"] = 2
    with pytest.raises(ValueError, match="outside"):
        aligned_reference(renderer, 1)
