"""An isolated attack spectrum cannot include energy from a later region."""

import importlib.util
from pathlib import Path
import numpy as np

spec = importlib.util.spec_from_file_location(
    "fit_plots", Path(__file__).parents[2] / "tools/instrument_fit_plots.py"
)
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)


def test_later_impulse_does_not_leak_into_attack_spectrum():
    rate = 48000
    samples = np.zeros(rate)
    samples[round(0.125 * rate)] = 1
    _, attack = plots.region_spectrum(samples, rate, 0, 0.12)
    _, later = plots.region_spectrum(samples, rate, 0.12, 0.5)
    assert not attack.any()
    assert later.any()
