"""Energy-conserving ERB-rate aggregation for spectrotemporal comparisons."""

from dataclasses import dataclass

import numpy as np


def erb_rate(frequency_hz: np.ndarray | float) -> np.ndarray:
    frequency = np.asarray(frequency_hz, dtype=np.float64)
    if np.any(frequency < 0) or not np.isfinite(frequency).all():
        raise ValueError("ERB frequencies must be finite and non-negative")
    return 21.4 * np.log10(1.0 + 0.00437 * frequency)


def frequency_from_erb_rate(rate: np.ndarray | float) -> np.ndarray:
    erb = np.asarray(rate, dtype=np.float64)
    if np.any(erb < 0) or not np.isfinite(erb).all():
        raise ValueError("ERB rates must be finite and non-negative")
    return (10.0 ** (erb / 21.4) - 1.0) / 0.00437


@dataclass(frozen=True)
class ErbFilterbank:
    frequencies_hz: np.ndarray
    centers_hz: np.ndarray
    weights: np.ndarray
    selected_bins: np.ndarray

    @classmethod
    def create(
        cls,
        frequencies_hz: np.ndarray,
        minimum_hz: float = 40.0,
        maximum_hz: float | None = None,
        band_count: int | None = None,
    ) -> "ErbFilterbank":
        frequencies = np.asarray(frequencies_hz, dtype=np.float64)
        if frequencies.ndim != 1 or frequencies.size < 2:
            raise ValueError("ERB filterbank needs a frequency-bin vector")
        maximum = float(frequencies[-1] if maximum_hz is None else maximum_hz)
        if not 0 <= minimum_hz < maximum <= frequencies[-1]:
            raise ValueError("ERB analysis range must lie inside the spectrum")
        minimum_rate = float(erb_rate(minimum_hz))
        maximum_rate = float(erb_rate(maximum))
        count = band_count or max(2, int(np.ceil(maximum_rate - minimum_rate)))
        if count < 2:
            raise ValueError("ERB filterbank needs at least two bands")
        center_rates = np.linspace(minimum_rate, maximum_rate, count)
        frequency_rates = erb_rate(frequencies)
        selected = (frequencies >= minimum_hz) & (frequencies <= maximum)
        weights = _partition_weights(frequency_rates, center_rates, selected)
        return cls(
            frequencies, frequency_from_erb_rate(center_rates), weights, selected
        )

    def apply_power(self, power: np.ndarray) -> np.ndarray:
        values = np.asarray(power, dtype=np.float64)
        if values.shape[0] != self.frequencies_hz.size:
            raise ValueError("spectral power does not match ERB frequency bins")
        if np.any(values < 0) or not np.isfinite(values).all():
            raise ValueError("spectral power must be finite and non-negative")
        return self.weights @ values


def _partition_weights(
    frequency_rates: np.ndarray, center_rates: np.ndarray, selected: np.ndarray
) -> np.ndarray:
    weights = np.zeros((center_rates.size, frequency_rates.size), dtype=np.float64)
    positions = np.searchsorted(center_rates, frequency_rates, side="right")
    for frequency_bin in np.flatnonzero(selected):
        right = int(positions[frequency_bin])
        if right == 0:
            weights[0, frequency_bin] = 1.0
        elif right >= center_rates.size:
            weights[-1, frequency_bin] = 1.0
        else:
            left = right - 1
            fraction = (frequency_rates[frequency_bin] - center_rates[left]) / (
                center_rates[right] - center_rates[left]
            )
            weights[left, frequency_bin] = 1.0 - fraction
            weights[right, frequency_bin] = fraction
    return weights
