"""Explicit unit-cube coordinates, constrained to the editable UI ranges."""

import numpy as np


class ParameterBox:
    def __init__(self, initial, bounds, descriptors):
        self.start, self.keys = dict(initial), list(bounds)
        if not self.keys:
            raise ValueError("At least one parameter is required")
        limits = np.array(list(bounds.values()), dtype=float)
        self.physical_limits = limits.copy()
        metadata = {d["key"]: d for d in descriptors}
        for key, (lo, hi) in zip(self.keys, limits):
            d = metadata[key]
            if not d["minimum"] <= lo < hi <= d["maximum"]:
                raise ValueError(f"Bounds outside UI: {key}")
            if not lo <= initial[key] <= hi:
                raise ValueError(f"Start outside shared experiment bounds: {key}")
        self.logarithmic = np.array(
            [metadata[k].get("scale", "linear") == "logarithmic" for k in self.keys]
        )
        if np.any(limits[self.logarithmic] <= 0):
            raise ValueError("Logarithmic bounds must be positive")
        limits[self.logarithmic] = np.log(limits[self.logarithmic])
        self.low, self.high = limits.T
        values = np.array([initial[k] for k in self.keys], dtype=float)
        values[self.logarithmic] = np.log(values[self.logarithmic])
        self.initial = (values - self.low) / (self.high - self.low)

    def unpack(self, coordinates):
        values = self.low + coordinates * (self.high - self.low)
        values[self.logarithmic] = np.exp(values[self.logarithmic])
        # exp(log(bound)) can exceed a strict DSP/API boundary by one ULP.
        values = np.clip(values, self.physical_limits[:, 0], self.physical_limits[:, 1])
        return dict(self.start, **dict(zip(self.keys, values.tolist())))
