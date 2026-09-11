"""Reference-anchored regional power and rise loss for layered percussion.

No audition normalization or synthetic target: every cell comes from the real
reference at its saved gain. Shape-only scoring is a screening approximation;
final candidates must match absolute levels with actual observation parameters.

This deliberately pooled objective cannot certify pitch/ridge preservation.
See the regional spectrum audit and the documented gong counterexamples before
using an improved aggregate score as evidence of a better instrument fit.
"""

import numpy as np
from scipy.signal import stft


class LayeredBandLoss:
    units = "reference band-envelope RMS dB"
    edges = np.array(
        [80, 250, 500, 800, 1250, 2000, 3200, 5000, 7000, 9000, 11500, 15000]
    )
    times = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 1.25, 1.5, 2, 3, 4, 6]

    def __init__(self, reference, rate, audibility=False):
        reference = np.asarray(reference)
        if (
            not np.isfinite(rate)
            or rate < 32000
            or reference.ndim != 1
            or len(reference) != round(6 * rate)
            or not np.isfinite(reference).all()
        ):
            raise ValueError(
                "Layered analysis requires six seconds of finite mono audio at >=32 kHz"
            )
        self.rate, self.frames = rate, len(reference)
        self.regions = list(zip(self.times[:-1], self.times[1:]))
        self.active = np.ones(len(self.edges) - 1, dtype=bool)
        self.floor = 1e-12
        self.target = self.envelopes(reference)
        # Explicit reference-only analysis floor, not independently scaled audio.
        self.lower = self.target.max(axis=1, keepdims=True) - 45
        self.target = np.maximum(self.target, self.lower)
        self.weights = np.ones_like(self.target)
        if audibility:
            relative = self.target - self.target.max(axis=1, keepdims=True)
            self.weights = np.maximum(0.05, 10 ** (relative / 20))
            self.weights /= self.weights.mean(axis=1, keepdims=True)
        self.specification = dict(
            version="reference-layered-band-v1",
            edges_hz=self.edges.tolist(),
            times_seconds=self.times,
            fft_size=4096,
            hop_seconds=0.01,
            floor="reference band peak minus 45 dB",
            normalization=False,
            objective="absolute cell error plus half-weight within-band shape error",
            audibility_weighting=(
                "reference-only amplitude ratio to each band peak, minimum .05; mean-one per band"
                if audibility
                else "uniform cells"
            ),
        )

    def envelopes(self, audio):
        audio = np.asarray(audio)
        if audio.shape != (self.frames,) or not np.isfinite(audio).all():
            raise ValueError("Expected finite reference-length mono audio")
        f, t, z = stft(
            audio, self.rate, nperseg=4096, noverlap=4096 - round(0.01 * self.rate)
        )
        bands = np.array(
            [
                np.sum(abs(z[(f >= a) & (f < b)]) ** 2, axis=0)
                for a, b in zip(self.edges[:-1], self.edges[1:])
            ]
        )
        power = np.array(
            [bands[:, (t >= a) & (t < b)].mean(axis=1) for a, b in self.regions]
        ).T
        return 10 * np.log10(np.maximum(self.floor, power))

    def residual(self, db, shape_only=False):
        db = np.asarray(db)
        if db.shape[-2:] != self.target.shape or not np.isfinite(db).all():
            raise ValueError("Expected finite matching band/time cells")
        error = np.maximum(db, self.lower) - self.target
        mean = (error * self.weights).mean(axis=-1, keepdims=True)
        shape = (error - mean) * np.sqrt(self.weights)
        error = error * np.sqrt(self.weights)
        if shape_only:
            return shape
        return np.concatenate([error, 0.5 * shape], axis=-1)

    def score_db(self, db, shape_only=False):
        errors = self.residual(db, shape_only)
        return float(np.sqrt(np.mean(errors**2)))

    def score(self, audio):
        return self.score_db(self.envelopes(audio))

    def diagnostics(self, audio):
        db = self.envelopes(audio)
        return dict(
            score=self.score_db(db),
            shape=self.score_db(db, True),
            envelopes_db=db.tolist(),
        )

    def attribution(self, audio):
        """Expose additive objective costs; a lower total need not improve each cell.

        Costs sum to score squared. Bias and shape remain separate diagnostics:
        reducing an observation gain cannot repair a wrong envelope shape.
        """
        db = self.envelopes(audio)
        error = np.maximum(db, self.lower) - self.target
        bias = (error * self.weights).mean(axis=-1, keepdims=True)
        cost = self.weights * (error**2 + 0.25 * (error - bias) ** 2)
        cost /= 2 * error.size
        return dict(
            score=float(np.sqrt(cost.sum())),
            cell_cost=cost.tolist(),
            band_cost=cost.sum(axis=1).tolist(),
            region_cost=cost.sum(axis=0).tolist(),
            band_bias_db=bias[:, 0].tolist(),
            band_shape_rms_db=np.sqrt(
                ((error - bias) ** 2 * self.weights).mean(axis=1)
            ).tolist(),
        )
