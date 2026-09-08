"""Absolute band/time energy subproblem, not a whole-sound acceptance score."""

import numpy as np
from scipy.signal import butter, sosfilt


class RegionalEnergyLoss:
    """Fit build-up and damping with fixed reference-derived floors.

    Integrate causal band-filter power inside explicit non-overlapping regions.
    No candidate normalization or peak alignment is performed. Frequency detail
    and within-region timing still require independent spectrogram checks.
    """

    units = "dB regional energy error"
    influence_threshold = 0.01

    def __init__(self, reference, rate, bands, regions):
        self.frames = len(reference)
        self.rate = float(rate)
        self.bands, self.regions = tuple(bands), tuple(regions)
        if not np.isfinite(rate) or rate <= 0 or not self.bands or not self.regions:
            raise ValueError("Expected positive rate and nonempty bands/regions")
        self.filters = []
        for low, high in self.bands:
            if not 0 < low < high < rate / 2:
                raise ValueError("Bands must lie strictly inside Nyquist")
            self.filters.append(
                butter(2, (low, high), fs=rate, btype="bandpass", output="sos")
            )
        self.slices = []
        previous = 0
        for low, high in self.regions:
            if not np.isfinite((low, high)).all() or not 0 <= low < high:
                raise ValueError("Expected finite nonnegative region boundaries")
            start, end = round(low * rate), round(high * rate)
            if not previous <= start < end <= self.frames:
                raise ValueError(
                    "Regions must be ordered, disjoint and inside the audio"
                )
            self.slices.append(slice(start, end))
            previous = end
        power = self.power(reference)
        self.floor = max(float(power.max()) * 1e-7, 1e-20)
        self.target = self.db(power)
        self.specification = dict(
            version="regional-energy-v1",
            bands_hz=self.bands,
            regions_seconds=self.regions,
            filter_order=2,
            causal=True,
            reference_floor_db=-70,
            normalization=False,
        )

    def power(self, samples):
        samples = np.asarray(samples)
        if (
            samples.shape != (self.frames,)
            or np.iscomplexobj(samples)
            or not np.isfinite(samples).all()
        ):
            raise ValueError("Expected finite mono audio matching the reference")
        rows = []
        for sos in self.filters:
            power = sosfilt(sos, samples) ** 2
            rows.append([power[section].mean() for section in self.slices])
        return np.array(rows)

    def db(self, power):
        return 10 * np.log10(np.maximum(power, self.floor))

    def residual(self, samples, regions=None):
        # Search's legacy five-region selector does not map to these regions.
        error = self.db(self.power(samples)) - self.target
        return error.ravel() / np.sqrt(error.size)

    def diagnostics(self, samples):
        difference = self.db(self.power(samples)) - self.target
        return dict(
            rms_error_db=float(np.sqrt(np.mean(difference**2))),
            bands_hz=self.bands,
            regions_seconds=self.regions,
            difference_db=difference.tolist(),
        )
