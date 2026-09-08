"""Exact quadratic band energies for an affine, actual-render audio basis."""

import numpy as np
from scipy.signal import sosfilt


class RegionalEnergyBasis:
    """Cache filtered cross-products; no rerender or FFT per amplitude step.

    Retain cross terms: summing individual modal powers would be incorrect for
    coherent signals, especially during the attack. This is an analysis cache,
    not an approximation of the instrument's internal energy states.
    """

    def __init__(self, basis, loss):
        self.loss = loss
        if not basis.bases or np.shape(basis.amplitudes) == (0,):
            raise ValueError("Energy basis needs at least one seed and amplitude")
        self.matrices = []
        for baseline, columns in basis.bases.values():
            if np.shape(columns) != (len(basis.amplitudes), loss.frames):
                raise ValueError("Basis columns must match amplitudes and duration")
            if (
                np.shape(baseline) != (loss.frames,)
                or not np.isfinite(columns).all()
                or not np.isfinite(baseline).all()
            ):
                raise ValueError("Expected finite matching basis waveforms")
            intercept = baseline - basis.amplitudes @ columns
            audio = np.vstack((intercept, columns))
            rows = []
            for sos in loss.filters:
                filtered = sosfilt(sos, audio, axis=-1)
                for section in loss.slices:
                    values = filtered[:, section]
                    rows.append(values @ values.T / values.shape[1])
            self.matrices.append(np.array(rows))
        self.matrices = np.array(self.matrices)
        self.normalization = np.sqrt(np.prod(self.matrices.shape[:2]))

    def evaluate(self, amplitudes):
        if np.shape(amplitudes) != (self.matrices.shape[-1] - 1,):
            raise ValueError("Wrong number of observation amplitudes")
        vector = np.r_[1.0, amplitudes]
        product = self.matrices @ vector
        power = product @ vector
        safe = np.maximum(power, self.loss.floor)
        error = 10 * np.log10(safe) - self.loss.target.ravel()
        jacobian = (20 / np.log(10)) * product[..., 1:] / safe[..., None]
        jacobian *= (power > self.loss.floor)[..., None]
        residual = error.ravel() / self.normalization
        jacobian = jacobian.reshape(-1, len(amplitudes)) / self.normalization
        if not np.isfinite(residual).all() or not np.isfinite(jacobian).all():
            raise ValueError("Nonfinite regional energy basis evaluation")
        return residual, jacobian

    def validate(self, basis, amplitudes):
        actual = np.concatenate(
            [
                self.loss.residual(audio + (amplitudes - basis.amplitudes) @ columns)
                for audio, columns in basis.bases.values()
            ]
        ) / np.sqrt(len(basis.bases))
        predicted = self.evaluate(amplitudes)[0]
        error = float(np.max(np.abs(actual - predicted)))
        if error > 1e-6:
            raise ValueError(
                f"Quadratic energy differs from waveform measurement: {error}"
            )
        return error


class RegionalEnergyGuard:
    """Bound each band/time error relative to an explicit comparator waveform.

    Analytic quadratic derivatives make this cheap enough for constrained
    observation fitting. The allowance is in dB, never a playback correction.
    """

    def __init__(self, basis, loss, comparator, tolerance_db=0.4):
        if not np.isfinite(tolerance_db) or tolerance_db < 0:
            raise ValueError("Expected nonnegative finite dB allowance")
        if set(comparator) != set(basis.bases):
            raise ValueError("Comparator must cover exactly the basis seeds")
        self.energy = RegionalEnergyBasis(basis, loss)
        self.limits = np.concatenate(
            [
                np.abs(loss.db(loss.power(comparator[seed])) - loss.target).ravel()
                + tolerance_db
                for seed in basis.bases
            ]
        )
        self.specification = dict(
            version="regional-energy-guard-v1",
            tolerance_db=tolerance_db,
            measurement=loss.specification,
            per_seed=True,
        )

    def evaluate(self, amplitudes):
        residual, jacobian = self.energy.evaluate(amplitudes)
        residual *= self.energy.normalization
        jacobian *= self.energy.normalization
        # Two smooth linear inequalities avoid abs()'s corner at zero error.
        margins = np.r_[self.limits - residual, self.limits + residual]
        return margins, np.vstack((-jacobian, jacobian))

    def values(self, amplitudes):
        return self.evaluate(amplitudes)[0]
