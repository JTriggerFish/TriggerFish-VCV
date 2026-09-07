"""Early RMS constraints over a validated, actual-render observation basis."""

import numpy as np

EARLY_BINS = ((0, 0.001), (0.001, 0.003), (0.003, 0.01), (0.01, 0.03), (0.03, 0.1))


class ObservationEnergyGate:
    """Return nonnegative constraint margins and exact amplitude Jacobians.

    Constraints apply separately to each training seed, never to averaged audio.
    A reference-only power floor protects silence without normalizing playback.
    """

    def __init__(self, basis, reference, rate, tolerance_db=3.5, bins=EARLY_BINS):
        if not np.isfinite(tolerance_db) or tolerance_db <= 0:
            raise ValueError("Energy tolerance must be finite and positive")
        if not np.isfinite(rate) or rate <= 0:
            raise ValueError("Sample rate must be finite and positive")
        reference = self._finite_vector(reference, "reference")
        initial = self._finite_vector(basis.amplitudes, "basis amplitudes")
        if not basis.bases:
            raise ValueError("At least one rendered basis is required")
        if getattr(basis, "sample_rate", rate) != rate:
            raise ValueError("Basis and reference sample rates differ")
        for audio, columns in basis.bases.values():
            audio = self._finite_vector(audio, "basis audio")
            columns = np.asarray(columns)
            if audio.shape != reference.shape or columns.shape != (
                len(initial),
                len(reference),
            ):
                raise ValueError("Basis and reference shapes differ")
            if not np.isfinite(columns).all():
                raise ValueError("Expected finite basis columns")
        self.frames = len(reference)
        self.basis, self.tolerance_db = basis, tolerance_db
        self.bins = tuple(bins)
        if not self.bins:
            raise ValueError("At least one energy bin is required")
        self.slices = []
        for start, end in self.bins:
            if not np.isfinite((start, end)).all():
                raise ValueError("Energy bin boundaries must be finite")
            first, last = round(start * rate), round(end * rate)
            if not 0 <= first < last <= len(reference):
                raise ValueError("Energy bin must lie inside the reference")
            self.slices.append(slice(first, last))
        power = self._power(reference)
        floor = max(float(np.max(power)) * 1e-10, 1e-20)
        self.floor = floor
        self.target = [max(float(np.mean(power[s])), floor) for s in self.slices]
        if not np.isfinite(self.target).all():
            raise ValueError("Nonfinite reference energy targets")

    @staticmethod
    def _finite_vector(values, name):
        values = np.asarray(values, dtype=np.float64)
        if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
            raise ValueError(f"Expected nonempty finite mono {name}")
        return values

    @staticmethod
    def _power(values):
        with np.errstate(over="ignore", invalid="ignore"):
            power = np.square(values)
        if not np.isfinite(power).all():
            raise ValueError("Energy calculation overflowed")
        return power

    def evaluate(self, amplitudes):
        """Return margins >= 0 for feasible amplitudes, and d(margin)/d(amplitude)."""
        amplitudes = self._finite_vector(amplitudes, "amplitudes")
        if amplitudes.shape != self.basis.amplitudes.shape:
            raise ValueError("Amplitude shape does not match basis")
        margins, rows = [], []
        for audio, columns in self.basis.bases.values():
            for section, target in zip(self.slices, self.target):
                local = columns[:, section]
                signal = audio[section] + (amplitudes - self.basis.amplitudes) @ local
                power = float(np.mean(self._power(signal)))
                error = 10 * np.log10(max(power, self.floor) / target)
                derivative = np.zeros(len(amplitudes))
                if power > self.floor:
                    derivative = (
                        20 / np.log(10) * (local @ signal) / len(signal) / power
                    )
                margins.extend((self.tolerance_db - error, self.tolerance_db + error))
                rows.extend((-derivative, derivative))
        margins, rows = np.asarray(margins), np.asarray(rows)
        if not np.isfinite(margins).all() or not np.isfinite(rows).all():
            raise ValueError("Nonfinite energy constraints or Jacobian")
        return margins, rows

    def actual_errors(self, samples):
        """Measure a fresh full render using the identical fixed reference targets."""
        samples = self._finite_vector(samples, "render")
        if len(samples) != self.frames:
            raise ValueError("Render and reference lengths differ")
        power = self._power(samples)
        errors = np.asarray(
            [
                10 * np.log10(max(float(np.mean(power[s])), self.floor) / target)
                for s, target in zip(self.slices, self.target)
            ]
        )
        if not np.isfinite(errors).all():
            raise ValueError("Nonfinite render energy errors")
        return errors
