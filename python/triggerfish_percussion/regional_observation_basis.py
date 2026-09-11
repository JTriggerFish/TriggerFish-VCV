"""Exact regional Welch power for an already-validated affine audio basis.

Frequency-cell integration matches RegionalSpectrumAudit, including partial FFT
bins. Preparation is off-line; evaluation is a small quadratic form per cell.
"""

import numpy as np
from scipy.signal import stft


class RegionalObservationBasis:
    def __init__(self, basis, audit, seed):
        baseline, columns = basis.bases[seed]
        intercept = baseline - basis.amplitudes @ columns
        audio = np.vstack((intercept, columns))
        matrices = []
        for start, end in audit.regions:
            segment = audio[:, round(start * audit.rate) : round(end * audit.rate)]
            size = min(8192, segment.shape[1])
            f, _, z = stft(
                segment,
                audit.rate,
                nperseg=size,
                noverlap=size // 2,
                nfft=32768,
                detrend="constant",
                boundary=None,
                padded=False,
                scaling="psd",
                axis=-1,
            )
            df = f[1] - f[0]
            one_sided = np.full(len(f), 2.0)
            one_sided[[0, -1]] = 1
            cells = []
            for lo, hi in zip(audit.edges[:-1], audit.edges[1:]):
                overlap = np.maximum(
                    0, np.minimum(f + df / 2, hi) - np.maximum(f - df / 2, lo)
                )
                active = overlap > 0
                values = (
                    z[:, active] * np.sqrt((overlap * one_sided)[active])[None, :, None]
                )
                values = values.reshape(len(audio), -1)
                cells.append((values @ values.conj().T).real / z.shape[-1])
            matrices.append(cells)
        self.matrices = np.array(matrices)

    def power(self, amplitudes):
        v = np.r_[1.0, amplitudes]
        return np.maximum(self.matrices @ v @ v, 0)

    def validate(self, basis, audit, seed, amplitudes):
        baseline, columns = basis.bases[seed]
        audio = baseline + (amplitudes - basis.amplitudes) @ columns
        predicted = audit.db(self.power(amplitudes))
        actual = audit.db(audit.power(audio))
        error = float(np.max(abs(predicted - actual)))
        if not np.isfinite(error) or error > 0.02:
            raise ValueError(f"Regional power cache mismatch: {error} dB")
        return error
