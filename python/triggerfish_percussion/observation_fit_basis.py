"""Exact-render affine observation basis, validated before using it in a search.

Only positive painted observation amplitudes may vary. Frequencies, active mode
membership, excitation, damping, turbulence and radiation must remain fixed.
The C++ field state is independent of observation gain; no DSP is reimplemented.
"""

import numpy as np


class ObservationBasis:
    def __init__(self, renderer, parameters, keys, seconds, seeds):
        if not keys or any(not key.startswith("resolved_level_") for key in keys):
            raise ValueError("Only active metallic observation levels are supported")
        if any(parameters[key] <= -71.99 for key in keys):
            raise ValueError("A basis cannot change active modal membership")
        self.renderer = renderer
        self.sample_rate = renderer.sample_rate
        self.metadata, self.initial = renderer.metadata, dict(parameters)
        self.keys, self.seconds = tuple(keys), seconds
        self.amplitudes = 10 ** (np.array([parameters[k] for k in keys]) / 20)
        self.bases, self.validation = {}, []
        for seed in seeds:
            self.bases[seed] = self.build(seed)
            self.validate(seed)

    def build(self, seed):
        baseline = self.renderer.render(self.initial, self.seconds, seed)
        columns = []
        for key, amplitude in zip(self.keys, self.amplitudes):
            # Probe an appreciable absolute amplitude, including for quiet bars.
            # A relative +6dB at -60dB produces a poor subtraction SNR when that
            # column is later raised by tens of dB during optimization.
            probe_db = 0 if self.initial[key] < -6 else -12
            changed = dict(self.initial, **{key: probe_db})
            actual = self.renderer.render(changed, self.seconds, seed)
            difference = 10 ** (changed[key] / 20) - amplitude
            columns.append((actual - baseline) / difference)
        return baseline, np.array(columns)

    def render(self, parameters, seconds, seed=None):
        if seconds != self.seconds or seed not in self.bases:
            raise ValueError("Basis duration or seed differs from its preparation")
        if parameters.keys() != self.initial.keys() or any(
            parameters[k] != v for k, v in self.initial.items() if k not in self.keys
        ):
            raise ValueError("A non-observation parameter changed")
        if any(not -71.99 < parameters[k] <= 6 for k in self.keys):
            raise ValueError("Observation outside positive active-mode range")
        baseline, columns = self.bases[seed]
        coefficients = 10 ** (np.array([parameters[k] for k in self.keys]) / 20)
        return baseline + (coefficients - self.amplitudes) @ columns

    def validate(self, seed):
        for sign in (-1, 1, 0):
            probe = dict(self.initial)
            for i, key in enumerate(self.keys):
                probe[key] = (
                    float(np.clip(probe[key] + sign * (1 + i % 5), -60, 6))
                    if sign
                    else float(-3 - i % 5)
                )
            predicted = self.render(probe, self.seconds, seed)
            actual = self.renderer.render(probe, self.seconds, seed)
            if not np.isfinite(predicted).all() or not np.isfinite(actual).all():
                raise ValueError("Observation validation requires finite audio")
            relative = float(
                np.linalg.norm(predicted - actual) / max(np.linalg.norm(actual), 1e-12)
            )
            maximum = float(np.max(np.abs(predicted - actual)))
            if (
                not np.isfinite((relative, maximum)).all()
                or relative > 3e-4
                or maximum > 3e-5
            ):
                raise ValueError(
                    f"Observation is not affine enough: {relative=}, {maximum=}"
                )
            self.validation.append(
                dict(seed=seed, relative_rms=relative, maximum_absolute=maximum)
            )

    def request(self, **request):
        return self.renderer.request(**request)
