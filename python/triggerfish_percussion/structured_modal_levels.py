"""Two fitting coordinates for an entire painted series, not one per ridge.

The incoming shape is preserved up to a common gain and tilt in dB/octave.
These are offline search coordinates, baked into the visible saved bar levels;
they add no hidden runtime controls. They cannot repair a poor starting series.
"""

import numpy as np


class StructuredModalLevels:
    def __init__(self, parameters):
        self.keys = tuple(
            f"resolved_level_{i}"
            for i in range(32)
            if parameters[f"resolved_level_{i}"] > -71.99
        )
        if not self.keys:
            raise ValueError("A structured level fit needs an active modal series")
        self.initial = np.array([parameters[k] for k in self.keys], dtype=float)
        frequencies = np.array(
            [parameters[k.replace("level", "frequency")] for k in self.keys],
            dtype=float,
        )
        if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
            raise ValueError("Modal frequencies must be finite and positive")
        if not np.all(np.isfinite(self.initial)) or np.any(self.initial > 6):
            raise ValueError("Modal levels exceed the visible control range")
        self.matrix = np.column_stack(
            (np.ones(len(self.keys)), np.log2(frequencies / 1000))
        )

    def levels(self, controls):
        """Apply [gain dB, tilt dB/octave], without clipping individual bars."""
        controls = np.asarray(controls, dtype=float)
        if controls.shape != (2,) or not np.all(np.isfinite(controls)):
            raise ValueError("Expected finite shared gain and tilt")
        return self.initial + self.matrix @ controls
