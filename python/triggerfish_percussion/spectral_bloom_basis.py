"""Exact STFT cross-power cache for fast, differentiable observation fitting."""

import numpy as np
from scipy.signal import stft


class SpectralBloomBasis:
    def __init__(self, basis, loss):
        self.loss = loss
        self.matrices = []
        for baseline, columns in basis.bases.values():
            intercept = baseline - basis.amplitudes @ columns
            audio = np.vstack((intercept, columns))
            f, t, z = stft(
                audio,
                loss.rate,
                nperseg=4096,
                noverlap=4096 - round(0.01 * loss.rate),
                boundary="zeros",
                axis=-1,
            )
            bands = []
            for lo, hi in zip(loss.edges[:-1], loss.edges[1:]):
                spectrum = z[:, (f >= lo) & (f < hi)]
                cells = []
                for start, end in loss.regions:
                    selected = spectrum[:, :, (t >= start) & (t < end)]
                    values = selected.reshape(len(audio), -1)
                    cells.append((values @ values.conj().T).real / selected.shape[-1])
                bands.append(cells)
            self.matrices.append(np.array(bands)[loss.active])
        self.matrices = np.array(self.matrices)

    def pack(self, errors):
        # Leading dimensions: seed, band, time; optional final derivative axis.
        count = np.prod(errors.shape[:3])
        rise = errors[:, :, 2:9] - errors[:, :, :1]
        trailing = errors.shape[3:]
        return np.concatenate(
            (
                errors.reshape((-1,) + trailing) / np.sqrt(count),
                rise.reshape((-1,) + trailing) / np.sqrt(np.prod(rise.shape[:3])),
            )
        )

    def evaluate(self, amplitudes):
        vector = np.r_[1.0, amplitudes]
        product = self.matrices @ vector
        power = product @ vector
        safe = np.maximum(power, self.loss.floor)
        error = 10 * np.log10(safe) - self.loss.target[self.loss.active]
        derivative = (20 / np.log(10)) * product[..., 1:] / safe[..., None]
        derivative *= (power > self.loss.floor)[..., None]
        return self.pack(error), self.pack(derivative)

    def validate(self, basis, amplitudes):
        # Compare each seed separately; pack uses component-major ordering.
        errors = np.array(
            [
                (
                    self.loss.db(
                        self.loss.power(
                            audio + (amplitudes - basis.amplitudes) @ columns
                        )
                    )
                    - self.loss.target
                )[self.loss.active]
                for audio, columns in basis.bases.values()
            ]
        )
        error = float(np.max(abs(self.pack(errors) - self.evaluate(amplitudes)[0])))
        if error > 1e-6:
            raise ValueError(f"STFT cross-power cache mismatch: {error}")
        return error


class SpectralBloomGuard:
    """Permit spectral polishing without giving back the fitted rise/decay."""

    def __init__(self, basis, loss, tolerance_db=1.0, baseline_audio=None):
        if not np.isfinite(tolerance_db) or tolerance_db < 0:
            raise ValueError("Expected nonnegative finite guard tolerance")
        self.cache = SpectralBloomBasis(basis, loss)
        baseline = self.cache.evaluate(basis.amplitudes)[0]
        if baseline_audio is not None:
            # A reduced curve need not represent the incoming modal bars.
            # Protect the actual incoming audio, not its coarse projection.
            errors = np.array(
                [
                    (loss.db(loss.power(audio)) - loss.target)[loss.active]
                    for audio in baseline_audio
                ]
            )
            if (
                errors.shape != self.cache.matrices.shape[:3]
                or not np.isfinite(errors).all()
            ):
                raise ValueError("Guard baseline must match the seed/band/time grid")
            baseline = self.cache.pack(errors)
        seeds, bands, regions = self.cache.matrices.shape[:3]
        envelope_count, rise_count = seeds * bands * regions, seeds * bands * 7
        self.limits = (
            abs(baseline)
            + np.r_[
                np.full(envelope_count, tolerance_db / np.sqrt(envelope_count)),
                np.full(rise_count, tolerance_db / np.sqrt(rise_count)),
            ]
        )
        self.specification = dict(
            version="spectral-bloom-guard-v1",
            tolerance_db=tolerance_db,
            measurement=loss.specification,
        )

    def evaluate(self, amplitudes):
        value, jacobian = self.cache.evaluate(amplitudes)
        return np.r_[self.limits - value, self.limits + value], np.vstack(
            (-jacobian, jacobian)
        )

    def values(self, amplitudes):
        return self.evaluate(amplitudes)[0]
