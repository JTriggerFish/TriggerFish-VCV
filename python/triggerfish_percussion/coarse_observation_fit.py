"""Broad observation controls, not independent noisy modal-bar fits.

Positive amplitudes interpolate in log frequency. This makes the actual output
affine in the coordinates and preserves a smooth, low-dimensional shape.
Coordinates are fitting tools only: the result stores ordinary visible bars.
"""

import numpy as np
from scipy.optimize import least_squares

from .spectral_bloom_basis import SpectralBloomBasis


def interpolation_weights(frequencies, knots):
    frequencies, knots = np.asarray(frequencies), np.asarray(knots)
    if (
        frequencies.ndim != 1
        or knots.ndim != 1
        or len(knots) < 2
        or not np.isfinite(frequencies).all()
        or not np.isfinite(knots).all()
        or np.any(frequencies <= 0)
        or np.any(knots <= 0)
        or np.any(np.diff(knots) <= 0)
    ):
        raise ValueError("Expected positive frequencies and ordered distinct knots")
    return np.array(
        [
            np.interp(np.log(frequencies), np.log(knots), row)
            for row in np.eye(len(knots))
        ]
    ).T


class CoarseObservationBasis:
    def __init__(self, renderer, parameters, seconds, seeds, knots):
        self.initial, self.renderer = dict(parameters), renderer
        self.keys = [
            f"resolved_level_{i}"
            for i in range(32)
            if parameters[f"resolved_level_{i}"] > -71.99
        ]
        frequencies = [
            parameters[key.replace("level", "frequency")] for key in self.keys
        ]
        self.weights = interpolation_weights(frequencies, knots)
        incoming = 10 ** (np.array([parameters[key] for key in self.keys]) / 20)
        projected = (
            self.weights @ np.linalg.lstsq(self.weights, incoming, rcond=None)[0]
        )
        if np.linalg.norm(incoming - projected) > 1e-6 * max(
            np.linalg.norm(incoming), 1e-12
        ):
            raise ValueError(
                "Source bars must already follow the coarse observation curve"
            )
        self.amplitudes = np.full(len(knots), 0.1)
        self.bases, self.validation = {}, []
        for seed in seeds:
            baseline = renderer.render(self.parameters(self.amplitudes), seconds, seed)
            columns = []
            for i in range(len(knots)):
                probe = self.amplitudes.copy()
                probe[i] = 1
                columns.append(
                    (renderer.render(self.parameters(probe), seconds, seed) - baseline)
                    / 0.9
                )
            columns = np.array(columns)
            probe = np.linspace(0.13, 1.2, len(knots))
            actual = renderer.render(self.parameters(probe), seconds, seed)
            predicted = baseline + (probe - self.amplitudes) @ columns
            error = float(
                np.linalg.norm(actual - predicted) / max(1e-12, np.linalg.norm(actual))
            )
            if not np.isfinite(error) or error > 3e-4:
                raise ValueError(f"Coarse observation basis is not affine: {error}")
            self.validation.append(dict(seed=seed, relative_error=error))
            self.bases[seed] = baseline, columns

    def parameters(self, amplitudes):
        amplitudes = np.asarray(amplitudes)
        if (
            amplitudes.shape != self.amplitudes.shape
            or not np.isfinite(amplitudes).all()
            or np.any(amplitudes <= 0)
            or np.any(amplitudes > 10 ** (6 / 20) + 1e-12)
        ):
            raise ValueError("Invalid coarse observation amplitudes")
        levels = 20 * np.log10(self.weights @ amplitudes)
        return dict(self.initial, **dict(zip(self.keys, levels.tolist())))


def polish_coarse(search, knots=(120, 600, 3000, 15000)):
    basis = CoarseObservationBasis(
        search.renderer, search.parameters, search.seconds, search.seeds, knots
    )
    cache = SpectralBloomBasis(basis, search.loss)
    low, high = 10 ** (-45 / 20), 10 ** (6 / 20)
    result = least_squares(
        lambda x: cache.evaluate(x)[0],
        basis.amplitudes,
        jac=lambda x: cache.evaluate(x)[1],
        bounds=(low, high),
        max_nfev=150,
    )
    validation = cache.validate(basis, result.x)
    candidate = basis.parameters(result.x)
    before = float(np.linalg.norm(search.residual(search.parameters)))
    after = float(np.linalg.norm(search.residual(candidate)))
    # Every incoming candidate must itself use this coarse shape. Never retain
    # a better score from an independently edited historical bar pattern.
    if after < before:
        search.parameters = candidate
    search.history.append(
        dict(
            stage=f"{len(knots)} broad observation amplitudes",
            knots_hz=list(knots),
            knots_db=(20 * np.log10(result.x)).tolist(),
            before=before,
            after=after,
            selected=after < before,
            basis_validation=basis.validation,
            cache_error=validation,
        )
    )
    search.save()
